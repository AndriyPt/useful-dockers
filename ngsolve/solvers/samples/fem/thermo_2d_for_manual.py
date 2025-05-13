# solve the Laplace equation -div lambda grad u = 0
# with Dirichlet boundary condition u = 0 and u = 2
# and Neumann side 0 conditions

from ngsolve import *
from netgen.occ import *
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


BORDER_TOP = "top"
BORDER_BOTTOM = "bottom"
BORDER_SIDE = "side"

DOMAIN_TISSUE = "tissue"
DOMAIN_TUMOR = "tumor"

INCLUSION_SIZE = 0.2
INCLUSION_HALF_SIZE = INCLUSION_SIZE / 2.0
INCLUSION_CENTER_X = 0.5
INCLUSION_CENTER_Y = 0.5
INCLUSION_RADIUS = INCLUSION_SIZE / 2.0 * math.sqrt(2.0)

K_TISSUE = 0.19  # W/m/^C
K_MAX_TUMOR = 0.495  # W/m/^C

LARGE_MAXH = 0.05

ngsglobals.msg_level = 1

whole_area = Rectangle(1.0, 1.0).Face()
whole_area.edges.Min(X).name = BORDER_SIDE
whole_area.edges.Min(Y).name = BORDER_BOTTOM
whole_area.edges.Max(X).name = BORDER_SIDE
whole_area.edges.Max(Y).name = BORDER_TOP

tumor_cross_section = (
    MoveTo(INCLUSION_CENTER_X - INCLUSION_HALF_SIZE, INCLUSION_CENTER_Y - INCLUSION_HALF_SIZE)
    .Rectangle(INCLUSION_SIZE, INCLUSION_SIZE)
    .Face()
)
tumor_cross_section.faces.name = DOMAIN_TUMOR

tissue_cross_section = whole_area - tumor_cross_section
tissue_cross_section.faces.name = DOMAIN_TISSUE

shape = Glue([tissue_cross_section, tumor_cross_section])

geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

print("Boundaries: ", mesh.GetBoundaries())
print("Materials: ", mesh.GetMaterials())

thermal_conductivity = mesh.MaterialCF(
    {
        DOMAIN_TISSUE: K_TISSUE,
        DOMAIN_TUMOR: (K_MAX_TUMOR - K_TISSUE)
        * cos(0.5 * pi / INCLUSION_RADIUS**2 * ((x - INCLUSION_CENTER_X) ** 2 + (y - INCLUSION_CENTER_Y) ** 2))
        + K_TISSUE,
    },
    default=0,
)

Draw(thermal_conductivity, mesh, "Thermal Conductivity")

# H1-conforming finite element space
fes = H1(mesh, order=3, dirichlet=BORDER_TOP + "|" + BORDER_BOTTOM)
dirichlet_condition = GridFunction(fes)
dirichlet_condition.Set(2.0 * y, BND)

# define trial- and test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# the right hand side
f = LinearForm(fes)

# the bilinear-form
a = BilinearForm(fes, symmetric=True)
a += thermal_conductivity * grad(u) * grad(v) * dx

a.Assemble()
f.Assemble()

# the solution field
res = f.vec.CreateVector()
res.data = f.vec - a.mat * dirichlet_condition.vec

gfu = GridFunction(fes)
gfu.vec.data = dirichlet_condition.vec.data + a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * res

# plot the solution (netgen-gui only)
Draw(gfu)

print("Visualizing data...")

CHART_STEPS = 20

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

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
        z_data[x_index][y_index] = gfu(point[0], point[1])

surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0)

plt.show()

print("Done!")
