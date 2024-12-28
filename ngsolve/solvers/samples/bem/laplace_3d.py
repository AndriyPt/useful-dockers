from netgen.occ import *
from ngsolve import *
from ngbem import *
from ngsolve.krylovspace import CG

sp = Sphere((0, 0, 0), 1)
mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=0.3)).Curve(3)
Draw(mesh)

fesL2 = SurfaceL2(mesh, order=3, dual_mapping=True)
u, v = fesL2.TnT()

u0 = 1 / sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)
Mu0 = LinearForm(u0 * v.Trace() * ds(bonus_intorder=3)).Assemble()

V = SingleLayerPotentialOperator(space=fesL2, intorder=10, method="aca")

j = GridFunction(fesL2)
pre = BilinearForm(u * v * ds, diagonal=True).Assemble().mat.Inverse()
with TaskManager(pajetrace=1000 * 1000 * 1000):
    CG(
        mat=V.mat,
        pre=pre,
        rhs=Mu0.vec,
        sol=j.vec,
        tol=1e-8,
        maxsteps=200,
        initialize=False,
        printrates=False,
    )
Draw(j)
