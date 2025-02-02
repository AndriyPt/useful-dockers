from netgen.occ import *
from ngsolve import *
from ngsolve import VTKOutput

sp = Sphere((0, 0, 0), 1)
mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=0.2)).Curve(4)

fes = H1(mesh, order=4, dirichlet=[1, 2, 3, 4])
u, v = fes.TnT()

bound = 1 / sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)

dirichlet_gfu = GridFunction(fes)
dirichlet_gfu.Set(bound, BND)

# the right hand side
f = LinearForm(fes)
f += 2.0 * v * dx

# the bilinear-form
a = BilinearForm(fes, symmetric=True)
a += grad(u) * grad(v) * dx

a.Assemble()
f.Assemble()

# the solution field
gfu = GridFunction(fes)
r = f.vec - a.mat * dirichlet_gfu.vec
gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * r + dirichlet_gfu.vec.data

vtk = VTKOutput(
    ma=mesh,
    coefs=[gfu, dirichlet_gfu],
    names=["temperature", "boundary"],
    filename="/tmp/thermo_fem_poisson_3d_result",
    subdivision=3,
)
vtk.Do()

# plot the solution (netgen-gui only)
Draw(gfu)
