from netgen.geom2d import *
from ngsolve import *

# This example is solving nonhomogeneous laplace equation
# with nonhomogeneous Dirichlet boundary conditions in 2D

geo = CSG2d()

size = 2.0
radius = 0.4

large_maxh = 0.5
small_maxh = 0.05

rect = Solid2d( [
    (0, 0),
    EdgeInfo(bc="right"),
    (size, 0),
    EdgeInfo(bc="top"),
    (size, size),
    EdgeInfo(bc="left"),
    (0, size),
    EdgeInfo(bc="bottom"),
  ], mat="outer_material" )

circle = Solid2d( [
    (0, -1),
    EdgeInfo(( 1,  -1), maxh=small_maxh), # control point for quadratic spline
    (1,0),
    EdgeInfo(( 1,  1), maxh=small_maxh), 
    (0,1),
    EdgeInfo((-1,  1), maxh=small_maxh),
    (-1,0),
    EdgeInfo((-1, -1), maxh=small_maxh),
    ], mat="inner_material")

circle.Scale(radius).Move((size / 2.0, size / 2.0))

domain1 = rect - circle
domain1.Mat("body_mat")
domain2 = circle
domain2.Mat("inner_mat")

geo.Add(domain1)
geo.Add(domain2)

mesh = Mesh(geo.GenerateMesh(maxh = large_maxh))
mesh.Curve(3)

# Draw(mesh)

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
