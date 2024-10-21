# Example from here https://ngsolve.org/forum/ngspy-forum/216-sontact-problems-of-the-theory-of-elasticity?start=0

from netgen.csg import *
from ngsolve import *

cube = OrthoBrick(Pnt(-1, -1, -1.5), Pnt(1,1,1)).bc("brick") \
  * Plane(Pnt(0,0,-1), Vec(0,0,-1)).bc("support")
cyl = Cylinder (Pnt(0,0,-1), Pnt(0,0,2), 0.3).bc("cyl") \
  * Plane(Pnt(0,0,0), Vec(0,0,-1)).bc("contact") \
  * Plane(Pnt(0,0,2), Vec(0,0,2)).bc("top")

geo = CSGeometry()
geo.Add( (cube-cyl).mat("mat1"), bcmod=[(cyl, "contact")] )
geo.Add ( cyl.mat("mat2"))

mesh = Mesh(geo.GenerateMesh(maxh=0.3))
mesh.Curve(2)

# Draw (mesh)
print ("bnds: ", mesh.GetBoundaries())
print ("mats: ", mesh.GetMaterials())

V1 = VectorH1(mesh, order=2, definedon="mat1", dirichlet="support")
V2 = VectorH1(mesh, order=2, definedon="mat2")
V = FESpace([V1,V2])

u1,u2 = V.TrialFunction()

g = CoefficientFunction((1,0,0))
a = BilinearForm(V)
a += SymbolicEnergy (0.5*InnerProduct(Sym(grad(u1)), Sym(grad(u1))), definedon=mesh.Materials("mat1"))
a += SymbolicEnergy (0.5*InnerProduct(Sym(grad(u2)), Sym(grad(u2))), definedon=mesh.Materials("mat2"))
a += SymbolicEnergy (g * u2, definedon=mesh.Boundaries("top"))


n = specialcf.normal(3)
def penalty(x):
    return IfPos(x, x*x, 0)
a += SymbolicEnergy ( penalty((u1-u2)*n), definedon=mesh.Boundaries("contact"))

gfu = GridFunction(V)
solvers.NewtonMinimization(a, gfu, inverse="sparsecholesky")

gfu1 = gfu.components[0]
gfu2 = gfu.components[1]
combinedu = CoefficientFunction( [gfu1, gfu2] )
# Draw (combinedu, mesh, "deformation")
# Draw ((gfu1-gfu2)*n, mesh, "gap")


# def penaltyprime(x):
#     return IfPos(x, 2*x, 0)
# contactforce = penaltyprime((gfu1-gfu2)*n)
# contactforce = CoefficientFunction( [contactforce if bc=="contact" else None for bc in mesh.GetBoundaries()] )
# Draw (contactforce, mesh, "contactforce")
