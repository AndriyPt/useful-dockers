from ngsolve import *
from netgen.geom2d import *

geo = CSG2d()

circle = Circle( center=(0,0), radius=1.0, mat="mat1", bc="bc_circle" )
rect = Rectangle( pmin=(0,0), pmax=(1.5,1.5), mat="mat2", bc="bc_rect" )

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
Draw(cf, mesh)


####################################


# geo = CSG2d()

# # define some primitives
# circle = Circle( center=(0,0), radius=1.0, mat="mat1", bc="bc_circle" )
# rect = Rectangle( pmin=(0,0), pmax=(1.5,1.5), mat="mat2", bc="bc_rect" )

# # use operators +, - and * for union, difference and intersection operations
# domain1 = circle - rect
# domain2 = circle * rect
# domain2.Mat("mat3").Maxh(0.1) # change domain name and maxh
# domain3 = rect-circle

# # add top level objects to geometry
# geo.Add(domain1)
# geo.Add(domain2)
# geo.Add(domain3)

# # generate mesh
# m = geo.GenerateMesh(maxh=0.3)

# mesh = Mesh(m)
# mesh.Curve(3)
# cf = mesh.RegionCF(VOL, dict(mat1=0, mat2=4, mat3=7))
# Draw(cf, mesh)

######################################

# geo = SplineGeometry()
# geo.AddRectangle ( (0,0), (2,2), bcs=["b","r","t","l"], leftdomain=1)
# geo.AddRectangle ( (1,0.9), (1.3,1.4), bcs=["b2","r2","t2","l2"], leftdomain=2, rightdomain=1)
# geo.SetMaterial (1, "outer")
# geo.SetMaterial (2, "inner")
# mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# fes1 = H1(mesh, definedon="inner")
# u1 = GridFunction(fes1, "u1")
# u1.Set (x*y)

# fes = H1(mesh, order=3,dirichlet="b|l|r")
# u = fes.TrialFunction()
# v = fes.TestFunction()

# gfu = GridFunction(fes)

# f = LinearForm(fes)
# f += SymbolicLFI (u1*v, definedon=mesh.Materials("inner"))
# f.Assemble()

# a = BilinearForm(fes)
# a += SymbolicBFI(grad(u)*grad(v))
# a.Assemble()

# gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

# Draw (u1)
# Draw (gfu)

#######################################

# ngsglobals.msg_level = 1

# geo = geom2d.SplineGeometry()
# pnums = [ geo.AddPoint (x,y,maxh=0.01) for x,y in [(0,0), (1,0), (1,0.1), (0,0.1)] ]
# for p1,p2,bc in [(0,1,"bot"), (1,2,"right"), (2,3,"top"), (3,0,"left")]:
#      geo.Append(["line", pnums[p1], pnums[p2]], bc=bc)
# mesh = Mesh(geo.GenerateMesh(maxh=0.05))


# # generate a triangular mesh of mesh-size 0.2
# # mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

# # H1-conforming finite element space
# fes = H1(mesh, order=2, dirichlet="left|right|top|bottom")

# # define trial- and test-functions
# u = fes.TrialFunction()
# v = fes.TestFunction()

# # the right hand side
# f = LinearForm(fes)
# f += 2 * v * dx

# # boundary conditions
# # 1. f'(right)*n = b, Neumann boundary
# neumann_b = -1
# f += neumann_b * v * ds(definedon='right')

# # the bilinear-form 
# a = BilinearForm(fes, symmetric=True)
# a += grad(u)*grad(v)*dx - 10 * u * v * dx

# a.Assemble()
# f.Assemble()

# # the solution field 
# gfu = GridFunction(fes)
# gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
# # print (u.vec)


# # plot the solution (netgen-gui only)
# Draw (gfu)
# Draw (-grad(gfu), mesh, "Flux")
