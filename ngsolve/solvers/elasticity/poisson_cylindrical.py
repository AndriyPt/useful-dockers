# https://math.stackexchange.com/questions/2832540/weak-form-of-poisson-equation-in-3d-cylindrical-space-heat-conduction
# Equation - \nabla u = -6
# Analytical solution in cylindrical coordinates r**2 + z**2  

from ngsolve import *
from netgen.occ import *

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

fes = H1(mesh, order=3, dirichlet=BORDER_BOTTOM + "|" + BORDER_TOP + "|" + BORDER_SIDE)

u = fes.TrialFunction()
v = fes.TestFunction()

f = LinearForm(fes)
f += -6 * v * dx

# the bilinear-form 
a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx

a.Assemble()
f.Assemble()

# the solution field 
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
# print (u.vec)


# plot the solution (netgen-gui only)
Draw (gfu)
Draw (-grad(gfu), mesh, "Flux")

exact = 16*x*(1-x)*y*(1-y)
print ("L2-error:", sqrt (Integrate ( (gfu-exact)*(gfu-exact), mesh)))


# # ------------------------------ LOAD LIBRARIES -------------------------------
# from netgen.geom2d import SplineGeometry
# from ngsolve import *
# from ngsolve.solvers import PreconditionedRichardson as PreRic
# from math import pi
# import sys

# # -------------------------------- PARAMETERS ---------------------------------
# h_max = 0.25
# k = 1

# inverse = "umfpack"

# # ----------------------------------- DATA ------------------------------------
# u_ex = CoefficientFunction((x * (1 - x) * y * (1 - y),
#                             (1 - x**2) * y * (1 - y)))
# x_rhs = CoefficientFunction( (x * 2 * y * (1 - y) + x * 2 * x * (1 - x)
#                              - (1 - 2 * x) * y * (1 - y)
#                              + (1 - x) * y * (1 - y),
#                             x * 4 * y * (1 - y) + x * 2 * (1 - x**2)))


# # ----------------------------------- MESH ------------------------------------
# geo = SplineGeometry()
# geo.AddRectangle((0, 0), (1, 1), bcs = ("bottom", "right", "top", "rot_axis"))
# mesh = Mesh(geo.GenerateMesh(maxh = h_max))

# # --------------------------- FINITE ELEMENT SPACE ----------------------------
# Vx = H1(mesh, order = k, dirichlet = "bottom|right|top|rot_axis")
# Vy = H1(mesh, order = k, dirichlet = "bottom|right|top")
# V = FESpace([Vx, Vy])

# freedofs = V.FreeDofs()
# gfu = GridFunction(V)

# # ----------------------------- (BI)LINEAR FORMS ------------------------------
# (ux, uy), (vx, vy) = V.TnT()

# grad_u = CoefficientFunction((grad(ux), grad(uy)), dims = (2,2))
# grad_v = CoefficientFunction((grad(vx), grad(vy)), dims = (2,2))

# diffusion = x * InnerProduct(grad_u, grad_v)
# diffusion += ux * vx / x
# forcing = x_rhs * CoefficientFunction((vx, vy))


# # -------------------------------- INTEGRATORS --------------------------------
# ir = IntegrationRule(points = [(1/3,1/3)], weights = [1] )

# a = BilinearForm(V)
# a += SymbolicBFI(diffusion).SetIntegrationRule(TRIG, ir)

# f = LinearForm(V)
# f += SymbolicLFI(forcing).SetIntegrationRule(TRIG, ir)


# # ---------------------------------- PROBLEM ----------------------------------
# a.Assemble()
# f.Assemble()

# from math import isnan
# nandofs = []

# rows,cols,vals = a.mat.COO()
# for i, j, v in zip(rows, cols, vals):
#     if isnan(v):
#         nandofs.extend([i, j])

# test = GridFunction(V)
# test.vec[:] = 0
# for i in nandofs:
#     test.vec[i] = 1

# Draw(CoefficientFunction((test.components[0],test.components[1])), mesh, "test")
