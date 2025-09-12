# solve the Pennes equation nabla u - k^2 u = 0
# with Dirichlet boundary condition u = 0 and u = 2
# analytical solution u = sinh(x)

from ngsolve import *
from netgen.occ import *

BORDER_TOP = "top"
BORDER_BOTTOM = "bottom"
BORDER_SIDE = "side"
LARGE_MAXH = 0.2
K_SQUARE = 1.0

ngsglobals.msg_level = 1

whole_area = Rectangle(1.0, 1.0).Face()
whole_area.edges.Min(X).name = BORDER_SIDE
whole_area.edges.Min(Y).name = BORDER_BOTTOM
whole_area.edges.Max(X).name = BORDER_SIDE
whole_area.edges.Max(Y).name = BORDER_TOP

shape = Glue([whole_area])

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)


# H1-conforming finite element space
fes = H1(mesh, order=3, dirichlet=BORDER_TOP + "|" + BORDER_BOTTOM + "|" + BORDER_SIDE)
dirichlet_condition = GridFunction(fes)
dirichlet_condition.Set(sinh(y), BND)

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)

# the bilinear-form
a = BilinearForm(fes, symmetric=True)
a += grad(u) * grad(v) * dx
a += -1.0 * K_SQUARE * u * v * dx

a.Assemble()
f.Assemble()

# the solution field
res = f.vec.CreateVector()
res.data = f.vec - a.mat * dirichlet_condition.vec

gfu = GridFunction(fes)
gfu.vec.data = dirichlet_condition.vec.data + a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res

# plot the solution (netgen-gui only)
# Draw(gfu)

exact = sinh(y)

error = GridFunction(fes)
error.Set(gfu - exact)
Draw(error)

print("L2-error:", sqrt(Integrate((gfu - exact) * (gfu - exact), mesh)))
