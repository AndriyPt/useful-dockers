from ngsolve import *
from netgen.occ import *


ROBIN_H = 4.184 # W/m^2/^C
ROBIN_ROOM_TEMPR = 12 # ^C

K_TISSUE = 0.19 # W/m/^C
K_MAX_TUMOR = 0.495 # W/m/^C
K_CO2 = 0.0176 # W/m/^C

TISSUE_METABOLIC_HEAT = 418.4 # W/m^3
ARTERY_TEMPR = 36.8 # ^C  
BLOOD_DENSITY = 1000 # kg/m^3
BLOOD_HEAT_CAPACITY = 4181 # J/kg/^C
TISSUE_PERFUSION = 0.0005 # m^3 / s / m^3

THERMAL_W = BLOOD_DENSITY * BLOOD_HEAT_CAPACITY * TISSUE_PERFUSION

BORDER_BOTTOM = "gamma_bottom"
BORDER_ENV = "gamma_environment"
BORDER_SIDE = "gamma_sides"

DOMAIN_TISSUE = "tissue"
DOMAIN_TUMOR = "tumor"
DOMAIN_CO2 = "co2"

PARAM_A = 0.02 # height (m)
PARAM_B = 0.16 # CO2 area height (m)
PARAM_C = 0.1 # width (m)

TUMOR_RADIUS = 0.0025 # m
TUMOR_CENTER_X = PARAM_C / 2.0
TUMOR_CENTER_Y = PARAM_A - TUMOR_RADIUS

LARGE_MAXH=0.005

# generate a triangular mesh
whole_area = Rectangle(PARAM_C, PARAM_A + PARAM_B).Face()
whole_area.edges.Min(X).name = BORDER_SIDE
whole_area.edges.Min(Y).name = BORDER_BOTTOM
whole_area.edges.Max(X).name = BORDER_SIDE
whole_area.edges.Max(Y).name = BORDER_ENV

bottom_area = Rectangle(PARAM_C, PARAM_A).Face()
co2_cross_section = whole_area - bottom_area 
co2_cross_section.faces.name = DOMAIN_CO2

tumor_cross_section = Circle((TUMOR_CENTER_X, TUMOR_CENTER_Y), TUMOR_RADIUS).Face()
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
    DOMAIN_TUMOR: (K_MAX_TUMOR - K_TISSUE) * cos(0.5 * pi / TUMOR_RADIUS ** 2 * ((x - TUMOR_CENTER_X) ** 2 + 
                                                 (y - TUMOR_CENTER_Y) ** 2)) + K_TISSUE,
    DOMAIN_CO2: K_CO2,
    }, 
    default = 0)

Draw(thermal_conductivity, mesh, "Thermal Conductivity")

heat_source = mesh.MaterialCF({ 
    DOMAIN_TISSUE: THERMAL_W * ARTERY_TEMPR + TISSUE_METABOLIC_HEAT,
    DOMAIN_TUMOR: THERMAL_W * ARTERY_TEMPR + TISSUE_METABOLIC_HEAT, # Should be OK for small Tumor
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
fes = H1(mesh, order=3, dirichlet=BORDER_BOTTOM)

# Dirichlet conditions
dirichlet_conditions = mesh.BoundaryCF({BORDER_BOTTOM: ARTERY_TEMPR}, default = 0)
dirichlet_gfu = GridFunction(fes)
dirichlet_gfu.Set(dirichlet_conditions, BND)


# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)
f += heat_source * v * dx + ROBIN_H * ROBIN_ROOM_TEMPR * v * ds(BORDER_ENV)

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
a += (thermal_conductivity * grad(u) * grad(v) + u_coeficient * u * v) * dx + ROBIN_H * u * v * ds(BORDER_ENV) 

a.Assemble()
f.Assemble()

# Approach for nonhomogeneous Dirichlet boundary condition
# https://docu.ngsolve.org/release/i-tutorials/unit-1.3-dirichlet/dirichlet.html
r = f.vec.CreateVector()
r.data = f.vec - a.mat * dirichlet_gfu.vec

# the solution field 
gfu = GridFunction(fes)
gfu.vec.data = dirichlet_gfu.vec.data + a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r

Draw(gfu, mesh, "Heat")
Draw(-grad(gfu), mesh, "Flux")
