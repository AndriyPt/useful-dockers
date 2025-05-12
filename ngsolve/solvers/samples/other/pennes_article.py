from ngsolve import *
from netgen.occ import *

RADIUS = 0.045 # m

ROBIN_H = 4.184 # W/m^2/^C
ROBIN_T_e = 12 # ^C

THERMAL_K = 0.6272 # W/m/^C

THERMAL_Q_m = 418.4 # W/m^3
THERMAL_T_a = 36.8 # ^C  
DENSITY_B = 1000 # kg/m^3
THERMAL_C_b = 4181 # J/kg/^C
THERMAL_omega = 0.0005 # m^3 / s / m^3

THERMAL_W = DENSITY_B * THERMAL_C_b * THERMAL_omega

BORDER_GAMMA = "gamma"
DOMAIN_TISSUE = "tissue"

LARGE_MAXH=0.005

ngsglobals.msg_level = 1

# generate a triangular mesh
arm_cross_section = Circle((0,0), RADIUS).Face()
arm_cross_section.edges.name = BORDER_GAMMA
arm_cross_section.faces.name = DOMAIN_TISSUE

shape = Glue([arm_cross_section])

geo = OCCGeometry(shape, dim = 2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

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
a += (THERMAL_K * grad(u) * grad(v) + THERMAL_W * u * v) * dx + ROBIN_H * u * v * ds(BORDER_GAMMA) 

a.Assemble()
f.Assemble()

# the solution field 
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

# plot the solution (netgen-gui only)
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")
