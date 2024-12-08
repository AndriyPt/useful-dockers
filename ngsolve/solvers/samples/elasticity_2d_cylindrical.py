# https://math.stackexchange.com/questions/2832540/weak-form-of-poisson-equation-in-3d-cylindrical-space-heat-conduction
# Equation - \nabla u = -6
# Analytical solution in cylindrical coordinates r**2 + z**2  

# NOT COMPLETED

from ngsolve import *
from netgen.occ import *

E, nu = 210, 0.2
mu  = E / 2 / (1 + nu)
lam = E * nu / ((1 + nu) * (1 - 2 * nu))

HEIGH = 0.10
RADIUS = 0.05

LARGE_MAXH=0.001

BORDER_ROTATION_AXIS = "rotation_axis"
BORDER_BOTTOM = "bottom"
BORDER_SIDE = "side"
BORDER_TOP = "top"

# generate a triangular mesh
whole_area = Rectangle(RADIUS, HEIGH).Face()
whole_area.edges.Min(X).name = BORDER_ROTATION_AXIS
whole_area.edges.Min(Y).name = BORDER_BOTTOM
whole_area.edges.Max(X).name = BORDER_SIDE
whole_area.edges.Max(Y).name = BORDER_TOP

shape = Glue([whole_area])

geo = OCCGeometry(shape, dim = 2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

def Stress(strain):
    return 2 * mu * strain + lam * Trace(strain) * Id(2)   

fes = VectorH1(mesh, order=3, dirichlet=BORDER_BOTTOM + "|" + BORDER_TOP + "|" + BORDER_SIDE)

# Dirichlet conditions
dirichlet_conditions = mesh.BoundaryCF({
    BORDER_BOTTOM: x * x, 
    BORDER_TOP: x * x + HEIGH * HEIGH,
    BORDER_SIDE: RADIUS * RADIUS + y * y}, 
    default = 0.0)
dirichlet_gfu = GridFunction(fes)
dirichlet_gfu.Set(dirichlet_conditions, BND)


u = fes.TrialFunction()
v = fes.TestFunction()

f = LinearForm(fes)
f += 2 * pi * (-6) * v * x * dx

a = BilinearForm(fes, symmetric=False)
a += 2 * pi * grad(u) * grad(v) * x * dx

# n = specialcf.normal(mesh.dim)
# a -= grad(u) * n * v * x * ds(definedon=BORDER_ROTATION_AXIS) # Should be zero since r = 0 on rotation axis

a.Assemble()
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * dirichlet_gfu.vec

gfu = GridFunction(fes)
gfu.vec.data = dirichlet_gfu.vec + a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r

Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")

exact = x*x + y*y
print ("L2-error:", sqrt(Integrate((gfu - exact) * (gfu - exact), mesh)))
