from ngsolve import *
from netgen.occ import *
from ngsolve import VTKOutput

# Bladder with tumor

ROBIN_H = 4.184 # W/m^2/^C
ROBIN_T_e = 21 # ^C m, Typical OR temperature

K_TISSUE = 0.19 # W/m/^C
K_MAX_TUMOR = 0.495 # W/m/^C
K_CO2 = 0.0176 # W/m/^C

THERMAL_TISSUE_Q_m = 1.21 * 1086 # W/m^3
THERMAL_TUMOR_Q_m = THERMAL_TISSUE_Q_m * 10 # W/m^3
THERMAL_T_a = 36.8 # ^C  
DENSITY_B = 1000 # kg/m^3
THERMAL_C_b = 4181 # J/kg/^C
THERMAL_TUMOR_omega = 0.0063 # m^3 / s / m^3
THERMAL_TISSUE_omega = 0.0031 # m^3 / s / m^3 

THERMAL_W = DENSITY_B * THERMAL_C_b * THERMAL_omega

BORDER_GAMMA = "gamma"
DOMAIN_TISSUE = "tissue"
DOMAIN_TUMOR = "tumor"
DOMAIN_CO2 = "co2"

PARAM_A = 0.04 # height (m)
PARAM_B = 0.16 # CO2 area height (m)
PARAM_C = 0.4 # width (m)

OMEGA3_RADIUS = 0.0025 # m
OMEGA3_CENTER_X = PARAM_C / 2.0
OMEGA3_CENTER_Y = PARAM_A - OMEGA3_RADIUS

LARGE_MAXH=0.005

# generate a triangular mesh
whole_area = Rectangle(PARAM_C, PARAM_A + PARAM_B).Face()
whole_area.edges.name = BORDER_GAMMA

bottom_area = Rectangle(PARAM_C, PARAM_A).Face()
co2_cross_section = whole_area - bottom_area 
co2_cross_section.faces.name = DOMAIN_CO2

tumor_cross_section = Circle((OMEGA3_CENTER_X, OMEGA3_CENTER_Y), OMEGA3_RADIUS).Face()
tumor_cross_section.faces.name = DOMAIN_TUMOR

tissue_cross_section = whole_area - co2_cross_section - tumor_cross_section 
tissue_cross_section.faces.name = DOMAIN_TISSUE

shape = Glue([tissue_cross_section, tumor_cross_section, co2_cross_section])

geo = OCCGeometry(shape, dim = 2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

thermal_conductivity = mesh.MaterialCF({
    DOMAIN_TISSUE: K_TISSUE,
    DOMAIN_TUMOR: K_MAX_TUMOR,
    DOMAIN_CO2: K_CO2,
    }, 
    default = 0)

Draw(thermal_conductivity, mesh, "Thermal Conductivity")

heat_source = mesh.MaterialCF({ 
    DOMAIN_TISSUE: THERMAL_W * THERMAL_T_a + THERMAL_Q_m,
    DOMAIN_TUMOR: THERMAL_W * THERMAL_T_a + THERMAL_Q_m * 10, # According to external article assumption
    },
    default = 0)

Draw(heat_source, mesh, "Heat Source")

u_coeficient = mesh.MaterialCF({ 
    DOMAIN_TISSUE: THERMAL_W,
    DOMAIN_TUMOR: THERMAL_W, # Should be OK for small Tumor
    },
    default = 0)

Draw(u_coeficient, mesh, "u Coefficient")


# H1-conforming finite element space
fes = H1(mesh, order=3)

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)
f += heat_source * v * dx + ROBIN_H * ROBIN_T_e * v * ds(BORDER_GAMMA)

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
a += (thermal_conductivity * grad(u) * grad(v) + u_coeficient * u * v) * dx + ROBIN_H * u * v * ds(BORDER_GAMMA) 

a.Assemble()
f.Assemble()

# the solution field 
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

# plot the solution (netgen-gui only)
Draw(gfu)

flux = -grad(gfu)
Draw(flux, mesh, "Flux")

# VTKOutput object
vtk = VTKOutput(ma=mesh,
                coefs=[gfu, flux],
                names = ["temperature", "flux"],
                filename="/tmp/thermo_result",
                subdivision=3)
vtk.Do()
