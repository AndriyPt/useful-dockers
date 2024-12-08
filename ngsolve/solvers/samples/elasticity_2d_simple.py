from netgen.occ import *
from ngsolve import *

LARGE_MAXH = 0.1

DIMENSIONS = 2

# generate a triangular mesh
whole_rect = Rectangle(3.0, 1.0).Face()
whole_rect.edges.name = "outer"
whole_rect.edges.Min(X).name = "fix"
whole_rect.edges.Max(X).name = "force"

holes = sum([Circle((0.5 + i, 0.5), 0.25).Face() for i in range(3)])
holes.edges.name = "cyl"

body = whole_rect - holes

shape = Glue([body])

geo = OCCGeometry(shape, dim = DIMENSIONS)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

E, nu = 210, 0.2
mu = E / 2 / (1 + nu)
lam = E * nu / ((1 + nu) * (1 - 2 * nu))

def Stress(strain):
    return 2 * mu * strain + lam * Trace(strain) * Id(DIMENSIONS)

fes = VectorH1(mesh, dim = DIMENSIONS, order = 3, dirichlet = "fix")
u, v = fes.TnT()
gfu = GridFunction(fes)

with TaskManager():
    a = BilinearForm(InnerProduct(Stress(Sym(Grad(u))), Sym(Grad(v))).Compile() * dx)
    pre = Preconditioner(a, "bddc")
    a.Assemble()

force = CF((1e-3, 0))
f = LinearForm(force * v * ds("force")).Assemble()

from ngsolve.krylovspace import CGSolver

inv = CGSolver(a.mat, pre, tol=1e-8)
gfu.vec.data = inv * f.vec

with TaskManager():
    fesstress = MatrixValued(H1(mesh, order=3), symmetric=True)
    gfstress = GridFunction(fesstress)
    gfstress.Interpolate(Stress(Sym(Grad(gfu))))

Draw(gfu, mesh, "displacement")

normal = CoefficientFunction((1.0, 0.0))

traction_x = gfstress * normal

Draw(traction_x, mesh, "traction_x_normal")
