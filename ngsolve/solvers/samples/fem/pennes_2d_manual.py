# solve the Pennes equation nabla u - k^2 u = 0
# with Dirichlet boundary condition u = 0 and u = 2
# analytical solution u = sinh(x)

from ngsolve import *
from netgen.geom2d import unit_square
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

LARGE_MAXH = 0.1
K_SQUARE = 1.0

ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size
mesh = Mesh(unit_square.GenerateMesh(maxh=LARGE_MAXH))

# H1-conforming finite element space
fes = H1(mesh, order=3, dirichlet=[1,2,3,4])

# exact = sinh(y)
exact = exp(y)

dirichlet_condition = GridFunction(fes)
dirichlet_condition.Set(exact, BND)

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)

# the bilinear-form
# a = BilinearForm(fes, symmetric=True)
a = BilinearForm(fes)
a += (grad(u) * grad(v) - K_SQUARE * u * v) * dx

a.Assemble()
f.Assemble()

# the solution field
res = f.vec.CreateVector()
res.data = f.vec - a.mat * dirichlet_condition.vec

gfu = GridFunction(fes)
gfu.vec.data = dirichlet_condition.vec.data + a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res

# plot the solution (netgen-gui only)
# Draw(gfu)

error = GridFunction(fes)
error.Set(gfu - exact)
Draw(error)

print("L2-error:", sqrt(Integrate((gfu - exact) * (gfu - exact), mesh)))

CHART_STEPS = 20

x_min = 0.0
y_min = 0.0
x_max = 1.0
y_max = 1.0

x_data_linear = np.linspace(x_min, x_max, CHART_STEPS)
y_data_linear = np.linspace(y_min, y_max, CHART_STEPS)

x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)

z_data = np.empty([CHART_STEPS, CHART_STEPS])
for x_index in range(CHART_STEPS):
    for y_index in range(CHART_STEPS):
        point = np.array([x_data_linear[x_index], y_data_linear[y_index]])
        z_data[x_index][y_index] = error(point[0], point[1])

if True:

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=10)

else:
    fig, ax = plt.subplots()

    CS = ax.contour(x_data, y_data, z_data, levels=10) # Draws 10 automatically chosen isolines

    ax.clabel(CS, inline=True, fontsize=8)

    # ax.set_title('Contour Plot of Z = sin(X) + cos(Y)')
    # ax.set_xlabel('X-axis')
    # ax.set_ylabel('Y-axis')

plt.show()

print("Done!")
