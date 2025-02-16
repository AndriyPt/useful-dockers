from netgen.occ import *
from ngsolve import *
from ngsolve import VTKOutput
from ngbem import *
from ngsolve.krylovspace import CG

POLINOM_ORDER = 3

sp = Sphere((0, 0, 0), 1)
sp.faces.name = "body"
sp.edges.name = "outer"
mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=0.2)).Curve(4)

# trace space, i.e. H^(1/2)(Gamma) conforming elements
fesH1 = H1(mesh, order=POLINOM_ORDER, definedon=mesh.Boundaries(".*"))
uH1, vH1 = fesH1.TnT()

# trace space, i.e. H^(-1/2)(Gamma) conforming elements
fesL2 = SurfaceL2(mesh, order=POLINOM_ORDER - 1, dual_mapping=True)
u, v = fesL2.TnT()
print("L2-ndof = ", fesL2.ndof, "LH1-ndof = ", fesH1.ndof)

uexa = 1 / sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)
u0 = GridFunction(fesH1)
u0.Interpolate(uexa)
Draw(u0, mesh, "origin")

Draw(uexa, mesh, "uexa")

# f = GridFunction(fesH1)
# f.Set(2.0)

# partial_sol = GridFunction(fesH1)
# u0_laplace = GridFunction(fesH1)
# with TaskManager():
#     V = SingleLayerPotentialOperator(fesL2, intorder=12, eps=1e-4)
#     partial_sol.Set(V.GetPotential(f), definedon=mesh.Boundaries(".*"), dual=False)
#     u0_laplace.Set(u0 - partial_sol)

# Draw(partial_sol, mesh, "partial_sol")

# Draw(u0_laplace, mesh, "u0_laplace")

# u1 = GridFunction(fesL2)
# pre = BilinearForm(u * v * ds, diagonal=True).Assemble().mat.Inverse()
# with TaskManager():
#     M = BilinearForm(uH1 * v.Trace() * ds(bonus_intorder=3)).Assemble()
#     K = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=12, leafsize=40, eta=3.0, eps=1e-4, method="aca")

#     rhs = ((0.5 * M.mat + K.mat) * u0_laplace.vec).Evaluate()
#     CG(
#         mat=V.mat,
#         pre=pre,
#         rhs=rhs,
#         sol=u1.vec,
#         tol=1e-8,
#         maxsteps=200,
#         initialize=False,
#         printrates=True,
#     )

# Draw(u1, mesh, "gradient")

# fesVolumeH1 = H1(mesh, order=POLINOM_ORDER)
# gf_screen = GridFunction(fesVolumeH1)
# print("ndofscreen=", fesVolumeH1.ndof)
# with TaskManager():
#     gf_screen.Set(
#         V.GetPotential(u1) - K.GetPotential(u0_laplace) + V.GetPotential(f),
#         definedon=mesh.Boundaries(".*"),
#         dual=False,
#     )
# Draw(gf_screen, mesh, "solution")

# vtk = VTKOutput(
#     ma=mesh,
#     coefs=[u1],
#     names=["temperature"],
#     filename="/tmp/thermo_bem_poisson_direct_3d_result",
#     subdivision=3,
# )
# vtk.Do()
