# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.solvers import PreconditionedRichardson as PreRic
from math import pi

# -------------------------------- PARAMETERS ---------------------------------
h_max = 0.25
k = 1

inverse = "umfpack"

# ----------------------------------- DATA ------------------------------------
u_ex = CoefficientFunction((x * (1 - x) * y * (1 - y),
                            (1 - x**2) * y * (1 - y)))
x_rhs = CoefficientFunction( (x * 2 * y * (1 - y) + x * 2 * x * (1 - x)
                             - (1 - 2 * x) * y * (1 - y)
                             + (1 - x) * y * (1 - y),
                            x * 4 * y * (1 - y) + x * 2 * (1 - x**2)))


# ----------------------------------- MESH ------------------------------------
geo = SplineGeometry()
geo.AddRectangle((0, 0), (1, 1), bcs = ("bottom", "right", "top", "rot_axis"))
mesh = Mesh(geo.GenerateMesh(maxh = h_max))

# --------------------------- FINITE ELEMENT SPACE ----------------------------
Vx = H1(mesh, order = k, dirichlet = "bottom|right|top|rot_axis")
Vy = H1(mesh, order = k, dirichlet = "bottom|right|top")
V = FESpace([Vx, Vy])

freedofs = V.FreeDofs()
gfu = GridFunction(V)

# ----------------------------- (BI)LINEAR FORMS ------------------------------
(ux, uy), (vx, vy) = V.TnT()

grad_u = CoefficientFunction((grad(ux), grad(uy)), dims = (2,2))
grad_v = CoefficientFunction((grad(vx), grad(vy)), dims = (2,2))

diffusion = x * InnerProduct(grad_u, grad_v)
diffusion += ux * vx / x
forcing = x_rhs * CoefficientFunction((vx, vy))


# -------------------------------- INTEGRATORS --------------------------------
ir = IntegrationRule(points = [(1/3,1/3)], weights = [1] )

a = BilinearForm(V)
a += SymbolicBFI(diffusion).SetIntegrationRule(TRIG, ir)

f = LinearForm(V)
f += SymbolicLFI(forcing).SetIntegrationRule(TRIG, ir)


# ---------------------------------- PROBLEM ----------------------------------
a.Assemble()
f.Assemble()

from math import isnan
nandofs = []

rows,cols,vals = a.mat.COO()
for i, j, v in zip(rows, cols, vals):
    if isnan(v):
        nandofs.extend([i, j])

test = GridFunction(V)
test.vec[:] = 0
for i in nandofs:
    test.vec[i] = 1

Draw(CoefficientFunction((test.components[0],test.components[1])), mesh, "test")
