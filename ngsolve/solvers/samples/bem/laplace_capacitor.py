# Source https://github.com/Weggler/ngbem/blob/main/demos/BEM_Laplace_Capacitor.ipynb

from ngsolve import *
from netgen.occ import *
import netgen.meshing as meshing
from ngsolve.krylovspace import CG, GMRes
from ngbem import *

largebox = Box((-2, -2, -2), (2, 2, 2))
b1 = Box((-1, -1, 0.5), (1, 1, 1))
b2 = Box((-1, -1, -1), (1, 1, -0.5))

largebox.faces.name = "outer"
b1.faces.name = "top"  # part of Gamma
b2.faces.name = "bot"  # part of Gamma
shell = largebox - b1 - b2  # Omega^c
shape = Compound([b1, b2])
mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=1))
Draw(mesh)

order = 3
fesH1 = H1(mesh, order=order, definedon=mesh.Boundaries(".*"))  # trace space, i.e. H^(1/2)(Gamma) conforming elements
uH1, vH1 = fesH1.TnT()
fesL2 = SurfaceL2(mesh, order=order - 1, dual_mapping=True)  # trace space, i.e. H^(-1/2)(Gamma) conforming elements
u, v = fesL2.TnT()
print("L2-ndof = ", fesL2.ndof, "LH1-ndof = ", fesH1.ndof)

utop = GridFunction(fesH1)
utop.Interpolate(1, definedon=mesh.Boundaries("top"))
ubot = GridFunction(fesH1)
ubot.Interpolate(-1, definedon=mesh.Boundaries("bot"))
u0 = utop.vec + ubot.vec

u1 = GridFunction(fesL2)
with TaskManager():
    V = SingleLayerPotentialOperator(fesL2, intorder=12, eps=1e-4)
    K = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=12, eps=1e-4)
    M = BilinearForm(uH1.Trace() * v.Trace() * ds(bonus_intorder=3)).Assemble()
    pre = BilinearForm(u.Trace() * v.Trace() * ds, diagonal=True).Assemble().mat.Inverse()
    rhs = ((-0.5 * M.mat + K.mat) * u0).Evaluate()
    CG(mat=V.mat, pre=pre, rhs=rhs, sol=u1.vec, tol=1e-8, maxsteps=200, initialize=False, printrates=False)

Draw(u1)

screen = (
    WorkPlane(Axes((0, 0, 0), X, Z)).RectangleC(4, 4).Face()
    - Box((-1.1, -1.1, 0.4), (1.1, 1.1, 1.1))
    - Box((-1.1, -1.1, -1.1), (1.1, 1.1, -0.4))
)
mesh_screen = Mesh(OCCGeometry(screen).GenerateMesh(maxh=0.25)).Curve(1)
fes_screen = H1(mesh_screen, order=3)
gf_screen = GridFunction(fes_screen)
print("ndofscreen=", fes_screen.ndof)
with TaskManager():
    gf_screen.Set(
        -V.GetPotential(u1) + K.GetPotential(utop) + K.GetPotential(ubot),
        definedon=mesh_screen.Boundaries(".*"),
        dual=False,
    )
Draw(gf_screen)

# Check with FEM

# mesh_shell = Mesh(OCCGeometry(shell).GenerateMesh(maxh=0.5)).Curve(1)
# fes_shell = H1(mesh_shell, order=3)
# gf_shell = GridFunction(fes_shell)
# with TaskManager():
#     gf_shell.Set(
#         -V.GetPotential(u1) + K.GetPotential(utop) + K.GetPotential(ubot),
#         definedon=mesh_shell.Boundaries("outer"),
#         dual=False,
#     )
# Draw(gf_shell)

# fesH1d = H1(mesh_shell, order=order, dirichlet="outer|top|bot")
# ud, vd = fesH1d.TnT()
# ad = BilinearForm(grad(ud) * grad(vd) * dx).Assemble()
# fd = LinearForm(fesH1d).Assemble()
# gfud = GridFunction(fesH1d)

# utop = GridFunction(fesH1d)
# utop.Interpolate(1, definedon=mesh_shell.Boundaries("top"))
# ubot = GridFunction(fesH1d)
# ubot.Interpolate(-1, definedon=mesh_shell.Boundaries("bot"))

# r = fd.vec - ad.mat * (gf_shell.vec + utop.vec + ubot.vec)
# gfud.vec.data = gf_shell.vec + utop.vec + ubot.vec
# gfud.vec.data += ad.mat.Inverse(freedofs=fesH1d.FreeDofs()) * r
#
# Draw(gfud)
