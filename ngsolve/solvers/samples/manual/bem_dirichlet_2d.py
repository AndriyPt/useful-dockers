#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

CHART_STEPS = 100
BORDER_ELEMENTS_COUNT = 10


class ExpressionTerm:
    pass


class BoundaryCondition:
    def __init__(self, border_elements: np.array):
        self.__border_elements = border_elements

    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def get_border_elements(self):
        return self.__border_elements


class BottomDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([0.0, 0.0]), np.array([1.0, 0.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 0.0


class LeftDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([0.0, 0.0]), np.array([0.0, 1.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 2.0 * point[1]


class RightDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([1.0, 0.0]), np.array([1.0, 1.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 2.0 * point[1]


class TopDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([0.0, 1.0]), np.array([1.0, 1.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 2.0


class ExpressionTerm:
    def __init__(self):
        pass

    def value(self, a: float, b: float):
        raise NotImplementedError("Call to abstract method")


class SingleLayerBoundaryTerm(ExpressionTerm):
    def __init__(self):
        super().__init__()

    def value(self, a: float, b: float):
        return a * b


class DoubleLayerBoundaryTerm(ExpressionTerm):
    def __init__(self):
        super().__init__()

    def value(self, a: float, b: float):
        return 10.0


class SingleLayerVolumeTerm(ExpressionTerm):
    def __init__(self):
        super().__init__()

    def value(self, a: float, b: float):
        return 10.0


print("Define boundary conditions...")

boundary_conditions = [
    BottomDirichletCondition(),
    LeftDirichletCondition(),
    RightDirichletCondition(),
    TopDirichletCondition(),
]

expression = [SingleLayerBoundaryTerm(), DoubleLayerBoundaryTerm(), SingleLayerVolumeTerm()]

unknown_count = 10

# print("Create SLAE...")

# matrix = None
# right_side = None

# for boundary in boundary_conditions:
#     for point in boundary.get_boundary_points():
#         matrix_row = np.zeros(unknown_count)
#         right_side_value = 0.0
#         for term in expression:
#             term.calculate_coefficients(matrix_row)
#             right_side += term.calculate_right_side()
#         if matrix is None:
#             matrix = matrix_row
#         else:
#             matrix = np.vstack((matrix, matrix_row))
#         if right_side is None:
#             right_side = np.array([right_side_value])
#         else:
#             right_side = np.append(right_side, [right_side_value])


# print("Solving SLAE...")

# solution = np.linalg.solve(matrix, right_side)

# print("Setting data back...")

# for term in expression:
#     term.propagate_solution(solution)

print("Visualizing data...")

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

x_min = min(np.min(boundary_condition.get_border_elements()[:, 0]) for boundary_condition in boundary_conditions)
y_min = min(np.min(boundary_condition.get_border_elements()[:, 1]) for boundary_condition in boundary_conditions)
x_max = max(np.max(boundary_condition.get_border_elements()[:, 0]) for boundary_condition in boundary_conditions)
y_max = max(np.max(boundary_condition.get_border_elements()[:, 1]) for boundary_condition in boundary_conditions)

x_data = np.linspace(x_min, x_max, CHART_STEPS)
y_data = np.linspace(y_min, y_max, CHART_STEPS)

x_data, y_data = np.meshgrid(x_data, y_data)

z_data = sum(np.vectorize(term.value)(x_data, y_data) for term in expression)

surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0)

plt.show()


print("Done!")
