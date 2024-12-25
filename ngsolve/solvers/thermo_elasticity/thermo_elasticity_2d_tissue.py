from netgen.occ import *
from ngsolve import *

E_TISSUE = 1060 # tissue
NU_TISSUE = 0.31

MU_TISSUE  = E_TISSUE / 2 / (1 + NU_TISSUE)
LAMBDA_TISSUE = E_TISSUE * NU_TISSUE / ((1 + NU_TISSUE) * (1 - 2 * NU_TISSUE))

DIMENSIONS = 2

BORDER_TOP = "gamma_top"
BORDER_PRESS = "gamma_pressure"
BORDER_SIDE = "gamma_side"
BORDER_BOTTOM = "gamma_bottom"

DOMAIN_TISSUE = "tissue"

PARAM_A = 0.04 # height (m)
PARAM_B = 0.16 # CO2 area height (m)
PARAM_C = 0.4 # width (m)

TOOL_RADIUS = 0.0025 # m
TOOL_DISPLACEMENT = -0.001 # m

OMEGA3_RADIUS = 0.0025 # m
OMEGA3_CENTER_X = PARAM_C / 2.0
OMEGA3_CENTER_Y = PARAM_A - OMEGA3_RADIUS

LARGE_MAXH=0.005

# generate a triangular mesh
left_rect = Rectangle(PARAM_C / 2.0 - TOOL_RADIUS, PARAM_A).Face()
left_rect.edges.Min(X).name = BORDER_SIDE
left_rect.edges.Max(Y).name = BORDER_TOP
left_rect.edges.Min(Y).name = BORDER_BOTTOM
left_rect.faces.name = DOMAIN_TISSUE

middle_rect = MoveTo(PARAM_C / 2.0 - TOOL_RADIUS, 0).Rectangle(2.0 * TOOL_RADIUS, PARAM_A).Face()
middle_rect.edges.Max(Y).name = BORDER_PRESS
middle_rect.edges.Min(Y).name = BORDER_BOTTOM
middle_rect.faces.name = DOMAIN_TISSUE

right_rect = MoveTo(PARAM_C / 2.0 + TOOL_RADIUS, 0).Rectangle(PARAM_C / 2.0 - TOOL_RADIUS, PARAM_A).Face()
right_rect.edges.Max(X).name = BORDER_SIDE
right_rect.edges.Max(Y).name = BORDER_TOP
right_rect.edges.Min(Y).name = BORDER_BOTTOM
right_rect.faces.name = DOMAIN_TISSUE

whole_body = Glue([left_rect, middle_rect, right_rect])
whole_body.faces.name = DOMAIN_TISSUE

shape = Glue([whole_body])

geo = OCCGeometry(shape, dim = DIMENSIONS)

mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

Draw(mesh)

lambda_coef = mesh.MaterialCF({
    DOMAIN_TISSUE: LAMBDA_TISSUE,
    }, 
    default = 0)

mu_coef = mesh.MaterialCF({
    DOMAIN_TISSUE: LAMBDA_TISSUE,
    }, 
    default = 0)


def Stress(strain):
    return 2 * mu_coef * strain + lambda_coef * Trace(strain) * Id(DIMENSIONS)   

fes = VectorH1(mesh, order=3, dirichlet=BORDER_BOTTOM + "|" + BORDER_SIDE + "|" + BORDER_PRESS)

# Dirichlet conditions
dirichlet_conditions = mesh.BoundaryCF({BORDER_PRESS: (0, TOOL_DISPLACEMENT)}, default = (0,0))
dirichlet_gfu = GridFunction(fes)
dirichlet_gfu.Set(dirichlet_conditions, BND)

u,v = fes.TnT()
gfu = GridFunction(fes)

with TaskManager():
    a = BilinearForm(InnerProduct(Stress(Sym(Grad(u))), Sym(Grad(v))).Compile()*dx)
    pre = Preconditioner(a, "bddc")
    a.Assemble()

from ngsolve.krylovspace import CGSolver
inv = CGSolver(a.mat, pre, tol=1e-8)

# Approach for nonhomogeneous Dirichlet boundary condition
# https://docu.ngsolve.org/release/i-tutorials/unit-1.3-dirichlet/dirichlet.html
f = LinearForm(fes) # Used to create vector of needed dimension 
r = f.vec.CreateVector()
r.data = inv * dirichlet_gfu.vec
gfu.vec.data = dirichlet_gfu.vec.data - inv * r

with TaskManager():
    fesstress = MatrixValued(H1(mesh, order=3), symmetric=True)
    gfstress = GridFunction(fesstress)
    gfstress.Interpolate(Stress(Sym(Grad(gfu))))

Draw(gfu, mesh, "Deformations")

normal = CoefficientFunction((0.0, 1.0))

traction = gfstress * normal 

Draw(traction, mesh, "Traction")

print ("Force applied (Newtons):", TOOL_RADIUS * Integrate(traction, mesh, definedon=mesh.Boundaries(BORDER_PRESS)))
