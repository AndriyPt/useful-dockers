from ngsolve import *
from netgen.occ import *


ROBIN_H = 4.184 # W/m^2/^C
ROBIN_T_e = 12 # ^C

THERMAL_K_tis = 0.19 # W/m/^C
THERMAL_K_tum = 0.495 # W/m/^C

THERMAL_Q_m = 418.4 # W/m^3
THERMAL_T_a = 36.8 # ^C  
DENSITY_B = 1000 # kg/m^3
THERMAL_C_b = 4181 # J/kg/^C
THERMAL_omega = 0.0005 # m^3 / s / m^3

THERMAL_W = DENSITY_B * THERMAL_C_b * THERMAL_omega

BORDER_GAMMA = "gamma"
DOMAIN_TISSUE = "tissue"
DOMAIN_TUMOR = "tumor"

PARAM_C = 0.4 # width (m)
PARAM_A = 0.04 # height (m)

OMEGA3_RADIUS = 0.0025 # m
OMEGA3_CENTER_X = PARAM_C / 2.0
OMEGA3_CENTER_Y = PARAM_A - OMEGA3_RADIUS

LARGE_MAXH=0.005

# generate a triangular mesh
whole_tissue = Rectangle(PARAM_C, PARAM_A).Face()
whole_tissue.edges.name = BORDER_GAMMA

tumor_cross_section = Circle((OMEGA3_CENTER_X, OMEGA3_CENTER_Y), OMEGA3_RADIUS).Face()
tumor_cross_section.faces.name = DOMAIN_TUMOR

tissue_cross_section = whole_tissue - tumor_cross_section 
tissue_cross_section.faces.name = DOMAIN_TISSUE

shape = Glue([tissue_cross_section, tumor_cross_section])

geo = OCCGeometry(shape, dim = 2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

thermal_conductivity = mesh.MaterialCF({
    DOMAIN_TISSUE: THERMAL_K_tis,
    DOMAIN_TUMOR: THERMAL_K_tum,
    }, 
    default = 0) 

# H1-conforming finite element space
fes = H1(mesh, order=5)

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)
f += (THERMAL_W * THERMAL_T_a + THERMAL_Q_m) * v * dx + ROBIN_H * ROBIN_T_e * v * ds(BORDER_GAMMA)

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
a += (thermal_conductivity * grad(u) * grad(v) + THERMAL_W * u * v) * dx + ROBIN_H * u * v * ds(BORDER_GAMMA) 

a.Assemble()
f.Assemble()

# the solution field 
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

# plot the solution (netgen-gui only)
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")
