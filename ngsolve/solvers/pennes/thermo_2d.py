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
BORDER_CONTACT = "contact"

DOMAIN_TISSUE = "tissue"
DOMAIN_CO2 = "co2"
DOMAIN_TUMOR = "tumor"

LARGE_MAXH = 0.1
SMALL_MAXH = 0.05

whole_box = Rectangle(PARAM_C, PARAM_A + PARAM_B).Face()
bottom_box = Rectangle(PARAM_C, PARAM_A).Face()

tumor_domain = Circle((OMEGA3_CENTER_X, OMEGA3_CENTER_Y), OMEGA3_RADIUS).Face()
tissue_domain = bottom_box - tumor_domain
co2_domain = whole_box - bottom_box  

tumor_domain.faces.name = DOMAIN_TUMOR
tissue_domain.faces.name = DOMAIN_TISSUE
co2_domain.faces.name = DOMAIN_CO2

tissue_domain.edges.Min(X).name = BORDER_GAMMA2
tissue_domain.edges.Min(Y).name = BORDER_GAMMA1
tissue_domain.edges.Max(X).name = BORDER_GAMMA2

co2_domain.edges.Max(Y).name = BORDER_GAMMA2
co2_domain.edges.Min(X).name = BORDER_GAMMA2
co2_domain.edges.Max(X).name = BORDER_GAMMA2

shape = Glue([tumor_domain, tissue_domain, co2_domain])

geo = OCCGeometry(shape, dim = 2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw (mesh)

fes = H1(mesh, order=3, dirichlet=BORDER_GAMMA2)
u = fes.TrialFunction()
v = fes.TestFunction()

gfu = GridFunction(fes)

thermal_conductivity = mesh.MaterialCF({
    DOMAIN_TISSUE: K_TISSUE,
    DOMAIN_TUMOR: (K_MAX_TUMOR - K_TISSUE) * cos(0.5 * pi / OMEGA3_RADIUS ** 2 * (x - OMEGA3_CENTER_X) ** 2 + 
                                                 (y - OMEGA3_CENTER_Y) ** 2) + K_TISSUE,
    DOMAIN_CO2: K_CO2
    }, 
    default = 0) 

heat_source = mesh.MaterialCF({ DOMAIN_TISSUE : 0.01 - (x - PARAM_C / 2) * (x - PARAM_C / 2) + (y - PARAM_A / 2) * (y - PARAM_A / 2) }, 
                              default = 0)

Draw(heat_source, mesh, "Source Function")

f = LinearForm(fes)
f += SymbolicLFI (heat_source * v)
f.Assemble()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u) * grad(v))
a.Assemble()

gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(gfu, mesh, "Heat")
