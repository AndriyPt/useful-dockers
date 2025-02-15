from netgen.occ import *
from ngsolve import *
from ngsolve import VTKOutput
from ngbem import *
from ngsolve.krylovspace import CG

sp = Sphere((0, 0, 0), 1)
sp.faces.name = "body"
sp.edges.name = "outer"
mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=0.2)).Curve(4)

fesL2 = SurfaceL2(mesh, order=3, dual_mapping=True)
u, v = fesL2.TnT()
fesH1 = H1(mesh, order=4)
uH1, vH1 = fesH1.TnT()
print("ndofL2 = ", fesL2.ndof, "ndof H1 = ", fesH1.ndof)

uexa = 1 / sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)
u0 = GridFunction(fesH1)
u0.Interpolate(uexa)
Draw(u0)

f = GridFunction(fesH1)
f.Set(2.0)

u1 = GridFunction(fesL2)
pre = BilinearForm(u * v * ds, diagonal=True).Assemble().mat.Inverse()
with TaskManager():
    V = SingleLayerPotentialOperator(
        fesL2,
        intorder=12,
        leafsize=40,
        eta=3.0,
        eps=1e-4,
        method="aca",
        testhmatrix=False,
    )

    M = BilinearForm(uH1 * v.Trace() * ds(bonus_intorder=3)).Assemble()
    K = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=12, leafsize=40, eta=3.0, eps=1e-4, method="aca")

    B = SingleLayerPotentialOperator(
        fesH1,
        definedon=mesh.Boundaries("outer"),
        intorder=12,
        leafsize=40,
        eta=3.0,
        eps=1e-4,
        method="aca",
        testhmatrix=False,
    )

    rhs = ((0.5 * M.mat + K.mat) * u0.vec + B.mat * f.vec).Evaluate()
    CG(
        mat=V.mat,
        pre=pre,
        rhs=rhs,
        sol=u1.vec,
        tol=1e-8,
        maxsteps=200,
        initialize=False,
        printrates=True,
    )

Draw(u1)

vtk = VTKOutput(
    ma=mesh,
    coefs=[u1],
    names=["temperature"],
    filename="/tmp/thermo_bem_poisson_direct_3d_result",
    subdivision=3,
)
vtk.Do()
