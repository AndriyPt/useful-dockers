from netgen.occ import *
from ngsolve import *

box = Box((0,0,0), (3,0.6,1))
box.faces.name="outer"
cyl = sum( [Cylinder((0.5+i,0,0.5), Y, 0.25,0.8) for i in range(3)] )
cyl.faces.name="cyl"
geo = box-cyl

cylboxedges = geo.faces["outer"].edges * geo.faces["cyl"].edges
cylboxedges.name = "cylbox"
geo = geo.MakeChamfer(cylboxedges, 0.03)

geo.faces.Min(X).name = "fix"
geo.faces.Max(X).name = "force"

mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.1)).Curve(3)
Draw(mesh)

E, nu = 210, 0.2
mu  = E / 2 / (1+nu)
lam = E * nu / ((1+nu)*(1-2*nu))

def Stress(strain):
    return 2*mu*strain + lam*Trace(strain)*Id(3)   

fes = VectorH1(mesh, order=3, dirichlet="fix")
u,v = fes.TnT()
gfu = GridFunction(fes)

with TaskManager():
    a = BilinearForm(InnerProduct(Stress(Sym(Grad(u))), Sym(Grad(v))).Compile()*dx)
    pre = Preconditioner(a, "bddc")
    a.Assemble()

force = CF( (1e-3,0,0) )
f = LinearForm(force*v*ds("force")).Assemble()

from ngsolve.krylovspace import CGSolver
inv = CGSolver(a.mat, pre, tol=1e-8)
gfu.vec.data = inv * f.vec

with TaskManager():
    fesstress = MatrixValued(H1(mesh,order=3), symmetric=True)
    gfstress = GridFunction(fesstress)
    gfstress.Interpolate (Stress(Sym(Grad(gfu))))

Draw(gfu, mesh, "Deformations")

Draw(Norm(gfstress), mesh, "Stress Norm")
