from netgen.geom2d import *
from ngsolve import *

a = 0.08
b = 0.16
c = 0.4

BORDER_GAMMA1 = "gamma1"
BORDER_GAMMA2 = "gamma2"
BORDER_CONTACT = "contact"

DOMAIN_TISSUE = "tissue"
DOMAIN_CO2 = "co2"

LARGE_MAXH = 0.1
SMALL_MAXH = 0.05

geo = SplineGeometry()
p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0, 0), (c, 0), (c, a), (0, a)] ]
p5,p6 =  [ geo.AppendPoint(x,y) for x,y in [(0, a + b), (c, a + b)] ]
geo.Append (["line", p1, p2], bc=BORDER_GAMMA1, leftdomain=1, rightdomain=0)
geo.Append (["line", p2, p3], bc=BORDER_GAMMA2, leftdomain=1, rightdomain=0)
geo.Append (["line", p3, p4], bc=BORDER_CONTACT, leftdomain=1, rightdomain=2)
geo.Append (["line", p4, p1], bc=BORDER_GAMMA2, leftdomain=1, rightdomain=0)
geo.Append (["line", p4, p5], bc=BORDER_GAMMA2, leftdomain=0, rightdomain=2)
geo.Append (["line", p5, p6], bc=BORDER_GAMMA2, leftdomain=0, rightdomain=2)
geo.Append (["line", p6, p3], bc=BORDER_GAMMA2, leftdomain=0, rightdomain=2)

geo.SetMaterial(1, DOMAIN_TISSUE)
geo.SetMaterial(2, DOMAIN_CO2)

mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH))

print ("bnds: ", mesh.GetBoundaries())
print ("mats: ", mesh.GetMaterials())

Draw(mesh)

fes = H1(mesh, order=3, dirichlet=BORDER_GAMMA2)
u = fes.TrialFunction()
v = fes.TestFunction()

gfu = GridFunction(fes)

heat_source = mesh.MaterialCF({ DOMAIN_TISSUE : 0.01 - (x - c / 2) * (x - c / 2) + (y - a / 2) * (y - a / 2) }, 
                              default = 0)

Draw(heat_source, mesh, "Source Function")

f = LinearForm(fes)
f += SymbolicLFI (heat_source * v)
f.Assemble()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u) * grad(v))
a.Assemble()

gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(gfu, mesh, "Heat")

sys.exit(0)

#################################3
V1 = VectorH1(mesh, order=2, definedon="mat1", dirichlet="support")
V2 = VectorH1(mesh, order=2, definedon="mat2")
V = FESpace([V1,V2])

u1,u2 = V.TrialFunction()

g = CoefficientFunction((1,0,0))
a = BilinearForm(V)
a += SymbolicEnergy (0.5*InnerProduct(Sym(grad(u1)), Sym(grad(u1))), definedon=mesh.Materials("mat1"))
a += SymbolicEnergy (0.5*InnerProduct(Sym(grad(u2)), Sym(grad(u2))), definedon=mesh.Materials("mat2"))
a += SymbolicEnergy (g * u2, definedon=mesh.Boundaries("top"))


fes1 = H1(mesh, definedon="inner")
u1 = GridFunction(fes1, "u1")
u1.Set (x*y)

fes = H1(mesh, order=3,dirichlet="b|l|r")
u = fes.TrialFunction()
v = fes.TestFunction()

gfu = GridFunction(fes)

f = LinearForm(fes)
f += SymbolicLFI (u1*v, definedon=mesh.Materials("inner"))
f.Assemble()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()

gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u1)
Draw (gfu)


sys.exit(0)
#################################

from netgen.geom2d import *
from ngsolve import *


geo = SplineGeometry()
geo.AddRectangle ( (0,0), (2,2), bcs=["b","r","t","l"], leftdomain=1)
geo.AddRectangle ( (1,0.9), (1.3,1.4), bcs=["b2","r2","t2","l2"], leftdomain=2, rightdomain=1)
geo.SetMaterial (1, "outer")
geo.SetMaterial (2, "inner")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

fes1 = H1(mesh, definedon="inner")
u1 = GridFunction(fes1, "u1")
u1.Set (x*y)

fes = H1(mesh, order=3,dirichlet="b|l|r")
u = fes.TrialFunction()
v = fes.TestFunction()

gfu = GridFunction(fes)

f = LinearForm(fes)
f += SymbolicLFI (u1*v, definedon=mesh.Materials("inner"))
f.Assemble()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()

gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u1)
Draw (gfu)

sys.exit(0)

#########################################

from ngsolve import *
from netgen.geom2d import *

# This example is solving homogeneous laplace equations
# on contacted rectangles from different materials in 2D

geo = CSG2d()

a = 0.08
b = 0.16
c = 0.4

large_maxh = 0.5
small_maxh = 0.05

tissue_rect = Solid2d( [
    (0, 0),
    EdgeInfo(bc="gamma2"),
    (c, 0),
    EdgeInfo(bc="contact"),
    (c, a),
    EdgeInfo(bc="gamma2"),
    (0, a),
    EdgeInfo(bc="gamma1"),
  ], mat="tissue_material" )

co2_rect = Solid2d( [
    (0, a),
    EdgeInfo(bc="gamma2"),
    (c, a),
    EdgeInfo(bc="gamma2"),
    (c, a + b),
    EdgeInfo(bc="gamma2"),
    (0, a + b),
    EdgeInfo(bc="contact"),
  ], mat="co2_material" )

geo.Add(tissue_rect)
geo.Add(co2_rect, bcmod=[(tissue_rect, "rect_contact")])

mesh = Mesh(geo.GenerateMesh(maxh = large_maxh))
mesh.Curve(3)

Draw(mesh)

sys.exit(0)

therm_conduct = cos(0.5 * pi / radius * ((x - 1.0) * (x - 1.0) + (y - 1.0) * (y - 1.0))) + 1.0
cf_therm_conduct = mesh.RegionCF(VOL, dict(body_mat=1, inner_mat=therm_conduct))

# Draw(cf_therm_conduct, mesh, "conductivity")

fes = H1(mesh, order=2, dirichlet="top|bottom")

# Dirichlet conditions
dirichlet_conditions = mesh.BoundaryCF({"bottom": 0, "top": 1}, default = 0)
gfu = GridFunction(fes)
gfu.Set(dirichlet_conditions, BND)

#Draw(gfu, mesh, "Dirichlet")

u = fes.TrialFunction()
v = fes.TestFunction()

f = LinearForm(fes)

a = BilinearForm(fes, symmetric=True)
a += cf_therm_conduct * grad(u) * grad(v) * dx

a.Assemble()
f.Assemble()

# Approach for nonhomogeneous Dirichlet boundary condition
# https://docu.ngsolve.org/release/i-tutorials/unit-1.3-dirichlet/dirichlet.html
r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec

# the solution field 
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

Draw (gfu, mesh, "Heat")
# Draw (-grad(gfu), mesh, "Flux")

sys.exit(0)

###################################

geo = CSG2d()

rect_co = Rectangle( pmin=(0,0), pmax=(1.5,1.5), mat="mat2", bc="bc_rect" )
rect_tissue = Rectangle( pmin=(0,0), pmax=(1.5,1.5), mat="mat2", bc="bc_rect" )

# use operators +, - and * for union, difference and intersection operations
domain1 = circle - rect
domain2 = circle * rect
domain2.Mat("mat3").Maxh(0.1) # change domain name and maxh
domain3 = rect-circle

# add top level objects to geometry
geo.Add(domain1)
geo.Add(domain2)
geo.Add(domain3)

# generate mesh
m = geo.GenerateMesh(maxh=0.3)

mesh = Mesh(m)
mesh.Curve(3)
cf = mesh.RegionCF(VOL, dict(mat1=0, mat2=4, mat3=7))

Draw(cf, mesh, "cf")


####################################

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
Draw (combinedu, mesh, "deformation")
Draw ((gfu1-gfu2)*n, mesh, "gap")


def penaltyprime(x):
    return IfPos(x, 2*x, 0)
contactforce = penaltyprime((gfu1-gfu2)*n)
contactforce = CoefficientFunction( [contactforce if bc=="contact" else None for bc in mesh.GetBoundaries()] )
Draw (contactforce, mesh, "contactforce")
