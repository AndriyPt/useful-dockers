from ngsolve import *
from netgen.occ import *


ROBIN_H = 4.184 # W/m^2/^C
ROBIN_T_e = 12 # ^C

K_TISSUE = 0.19 # W/m/^C
K_MAX_TUMOR = 0.495 # W/m/^C
K_CO2 = 0.0176 # W/m/^C

THERMAL_Q_m = 418.4 # W/m^3
THERMAL_T_a = 36.8 # ^C  
DENSITY_B = 1000 # kg/m^3
THERMAL_C_b = 4181 # J/kg/^C
THERMAL_omega = 0.0005 # m^3 / s / m^3

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
    DOMAIN_TUMOR: THERMAL_W * THERMAL_T_a + THERMAL_Q_m, # Should be OK for small Tumor
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
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")

sys.exit(0)

##########################################


from ngsolve import *
from netgen.occ import *

PARAM_A = 0.08 # m
PARAM_B = 0.16 # m
PARAM_C = 0.4 # m

OMEGA3_RADIUS = 0.005 # m
OMEGA3_CENTER_X = PARAM_C / 2.0
OMEGA3_CENTER_Y = PARAM_A - OMEGA3_RADIUS

K_TISSUE = 0.19 # W/m/K ==  W/m/Celcius
K_MAX_TUMOR = 0.495 # W/m/Celcius
K_CO2 = 0.0176 # W/m/Celcius

PERFUSION_OMEGA = 1.5 / 1000.0 / 60.0 * 10.0 # l/s/kg
BLOOD_DENSITY = 1050 # kg/m**3
BLOOD_HEAT_CAPACITY = 3617 # J/kg/Celcius
BODY_TEMP = 37 # Celcius
TISSUE_HEAT_GENERATION = 1.2 * 1086  # W/m**3, based on data Urinary Bladder Wall 1.21 W/kg, density 1086 kg/m**3

BORDER_GAMMA1 = "gamma1"
BORDER_GAMMA2 = "gamma2"
BORDER_GAMMA3 = "gamma3"

DOMAIN_TISSUE = "tissue"

ROBIN_BC_APLHA = 5
ROBIN_BC_ROOM_TEMP = 23

LARGE_MAXH = 0.1
SMALL_MAXH = 0.05

whole_box = Rectangle(PARAM_C, PARAM_A + PARAM_B).Face()

whole_box.faces.name = DOMAIN_TISSUE

whole_box.edges.Min(X).name = BORDER_GAMMA2
whole_box.edges.Min(Y).name = BORDER_GAMMA1
whole_box.edges.Max(X).name = BORDER_GAMMA2
whole_box.edges.Max(Y).name = BORDER_GAMMA3

shape = Glue([whole_box])

geo = OCCGeometry(shape, dim = 2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

fes = H1(mesh, order=3, dirichlet=BORDER_GAMMA1)
u = fes.TrialFunction()
v = fes.TestFunction()

# Dirichlet conditions
dirichlet_conditions = mesh.BoundaryCF({BORDER_GAMMA1: BODY_TEMP}, default = 0)
gfu = GridFunction(fes)
gfu.Set(dirichlet_conditions, BND)

Draw(gfu, mesh, "Dirichlet substitution")

thermal_conductivity = mesh.MaterialCF({
    DOMAIN_TISSUE: K_TISSUE,
    # DOMAIN_TUMOR: (K_MAX_TUMOR - K_TISSUE) * cos(0.5 * pi / OMEGA3_RADIUS ** 2 * ((x - OMEGA3_CENTER_X) ** 2 + 
    #                                              (y - OMEGA3_CENTER_Y) ** 2)) + K_TISSUE,
    # DOMAIN_CO2: K_CO2
    }, 
    default = 0) 

Draw(thermal_conductivity, mesh, "Thermal Conductivity")

heat_source = mesh.MaterialCF({ 
    DOMAIN_TISSUE: PERFUSION_OMEGA * BLOOD_DENSITY * BLOOD_HEAT_CAPACITY * BODY_TEMP + TISSUE_HEAT_GENERATION,
    # DOMAIN_TUMOR: PERFUSION_OMEGA * BLOOD_DENSITY * BLOOD_HEAT_CAPACITY * BODY_TEMP + TISSUE_HEAT_GENERATION  # Could be OK for small Tumor
    },
    default = 0)

Draw(heat_source, mesh, "Heat Source")

u_coeficient = mesh.MaterialCF({ 
    DOMAIN_TISSUE: -1.0 * PERFUSION_OMEGA * BLOOD_DENSITY * BLOOD_HEAT_CAPACITY,
    # DOMAIN_TUMOR: -1.0 * PERFUSION_OMEGA * BLOOD_DENSITY * BLOOD_HEAT_CAPACITY  # Could be OK for small Tumor
    },
    default = 0)

Draw(u_coeficient, mesh, "u Coefficient")

a = BilinearForm(fes, symmetric=True)
a += thermal_conductivity * grad(u) * grad(v) * dx + u_coeficient * u * v * dx + ROBIN_BC_APLHA * u * v * ds(BORDER_GAMMA3)

f = LinearForm(fes)
f += heat_source * v * dx + ROBIN_BC_APLHA * (ROBIN_BC_ROOM_TEMP - BODY_TEMP) * v * ds(BORDER_GAMMA3)

f.Assemble()
a.Assemble()

# Approach for nonhomogeneous Dirichlet boundary condition
# https://docu.ngsolve.org/release/i-tutorials/unit-1.3-dirichlet/dirichlet.html
r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec

# the solution field 
gfu_sol = GridFunction(fes)
gfu_sol.vec.data = gfu.vec.data + a.mat.Inverse(freedofs=fes.FreeDofs()) * r

Draw(gfu_sol, mesh, "Heat")
