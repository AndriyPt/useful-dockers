# solve the Pennes equation nabla u - k^2 u = 0
# with Dirichlet boundary condition u = 0 and u = 2
# analytical solution u = sinh(x)

from ngsolve import *
from netgen.geom2d import CSG2d, EdgeInfo as EI, PointInfo as PI, Solid2d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, path


LARGE_MAXH = 0.2
K_SQUARE = 0.25

ngsglobals.msg_level = 1

geo = CSG2d()

vertices = [
    (1, 0),
    (5, 0),
    (6, 1),
    (5, 2),
    (1, 2),
    (0, 1),
]

rect = Solid2d(vertices, mat="tissue")

geo.Add(rect)

# generate a triangular mesh of mesh-size
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH))

Draw(mesh)

# H1-conforming finite element space
fes = H1(mesh, order=3, dirichlet=[1, 2, 3, 4, 5, 6])

exact = sinh(0.5 * y)
# exact = exp(0.5 * y)

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

CHART_STEPS = 50

x_min = 0.0
y_min = 0.0
x_max = 6.0
y_max = 2.0

domain = path.Path(vertices)

x_data_linear = np.linspace(x_min, x_max, CHART_STEPS)
y_data_linear = np.linspace(y_min, y_max, CHART_STEPS)

x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)

z_data = np.empty([CHART_STEPS, CHART_STEPS])
z_data = np.zeros_like(x_data)
z_data.fill(np.nan)

for i in range(x_data.shape[0]):
    for j in range(x_data.shape[1]):
        point = (x_data[i, j], y_data[i, j])
        if domain.contains_point(point):
            z_data[i, j] = error(point[0], point[1])

if True:
    if True:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        fig.colorbar(surf, shrink=0.5, aspect=10)

    else:
        fig, ax = plt.subplots()

        CS = ax.contour(x_data, y_data, z_data, levels=10)  # Draws 10 automatically chosen isolines

        ax.clabel(CS, inline=True, fontsize=8)
    plt.show()

if False:
    import csv

    STEPS_CSV = 5
    x_values = np.linspace(0, 4, num=STEPS_CSV)
    y_values = np.linspace(0, 2, num=STEPS_CSV)

    header = ["/"]
    for y in y_values:
        header.append(y)

    data = []
    for x in x_values:
        line = [x]
        for y in y_values:
            value = abs(error(x, y))
            line.append(f"{value:.3}")
        data.append(line)

    # Specify the file path where the CSV will be saved
    file_path = "/home/user/workspace/project/ngsolve/solvers/samples/manual/fem_hexagon.csv"

    # Write data to CSV
    with open(file_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)  # Write the header first
        writer.writerows(data)  # Write the function values

    print(f"Data written to {file_path}")

print("Done!")
