from netgen.geom2d import *
from ngsolve import *

param_a = 0.08
param_b = 0.16
param_c = 0.4

BORDER_GAMMA1 = "gamma1"
BORDER_GAMMA2 = "gamma2"
BORDER_CONTACT = "contact"

DOMAIN_TISSUE = "tissue"
DOMAIN_CO2 = "co2"

LARGE_MAXH = 0.1
SMALL_MAXH = 0.05

geo = SplineGeometry()
p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0, 0), (param_c, 0), (param_c, param_a), (0, param_a)] ]
p5,p6 =  [ geo.AppendPoint(x,y) for x,y in [(0, param_a + param_b), (param_c, param_a + param_b)] ]
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

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

fes = H1(mesh, order=3, dirichlet=BORDER_GAMMA2)
u = fes.TrialFunction()
v = fes.TestFunction()

gfu = GridFunction(fes)

heat_source = mesh.MaterialCF({ DOMAIN_TISSUE : 0.01 - (x - param_c / 2) * (x - param_c / 2) + (y - param_a / 2) * (y - param_a / 2) }, 
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
