#!/usr/bin/env python3

from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

CHART_STEPS = 100
BORDER_ELEMENTS_COUNT = 10


class ExpressionTerm:
    pass

class BoundaryConditionType(Enum):
    DIRICHLET = 1
    NEUMANN = 2
    ROBIN = 3
    UNKNOWN = 4

class Point2DInfo:
    def __init__(self):
        self.point = np.zeros(2)
        self.type = BoundaryConditionType.UNKNOWN
        self.normal = np.zeros(2)
        self.element = np.empty(0)

class BoundaryCondition:
    def __init__(self, border_elements: np.array, type: BoundaryConditionType):
        assert border_elements is not None
        self.__border_elements = border_elements
        for point in self 


    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def get_border_elements(self):
        return self.__border_elements

    def get_boundary_points(self):
        if len(self.__border_elements) > 1:
            return self.__border_elements + (self.__border_elements[1] - self.__border_elements[0]) * 0.5
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


class Kernel:
    def value(self, point_x: np.array, point_y: np.array, normal: np.array):
        raise NotImplementedError("Call to abstract method")

class Laplace2DKernel:

    def value(self, point_x: np.array, point_y: np.array, normal: np.array):
        temp = point_x - point_y
        distance_sq = np.dot(temp.T, temp)
        return 0.25 / np.pi * np.log(distance_sq)

class Laplace2DNormKernel:

    def value(self, point_x: np.array, point_y: np.array, normal: np.array):
        raise NotImplementedError("Implement")
    

class Integrator:
    def calculate(self, kernel: Kernel, normal: np.array, point_x: np.array, point_y: np.array):
        # TODO: Implement
        return 0.0

class ExpressionTerm:
    def __init__(self, kernel: Kernel):
        self.__kernel = kernel
        assert self.__kernel is not None

    def value(self, a: float, b: float):
        raise NotImplementedError("Call to abstract method")
    
    def kernel(self):
        return self.__kernel

    def calculate_coefficients(
        self, point: np.array, boundary_conditions: list[BoundaryCondition], matrix_row: np.array
    ):
        raise NotImplementedError("Call to abstract method")

    def calculate_right_side(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        raise NotImplementedError("Call to abstract method")


class SingleLayerBoundaryTerm(ExpressionTerm):
    def __init__(self, kernel: Kernel):
        super().__init__(kernel)

    def calculate_coefficients(
        self, point: np.array, boundary_conditions: list[BoundaryCondition], matrix_row: np.array
    ):
        integrator = Integrator()
        coefficients = np.empty(0)
        right_side_ret = 0.0

        for boundary_condition in boundary_conditions:
            for boundary_point, element, normal, type in boundary_condition.get_boundary_points_info():
                if BoundaryConditionType.DIRICHLET == type:
                    res = integrator.calculate(self.kernel(), normal, point, boundary_point)
                    coefficients = np.append(coefficients, [res])
                elif BoundaryConditionType.NEUMANN == type:
                    right_side_ret += integrator.calculate(self.kernel(), normal, point, boundary_point)
                elif BoundaryConditionType.ROBIN == type:
                    raise NotImplementedError("Not implemented")
                else:
                    raise AttributeError("Not supported boundary element type")
        
        return (coefficients, right_side_ret)

    def calculate_right_side(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        pass

    def value(self, a: float, b: float):
        return a * b


class DoubleLayerBoundaryTerm(ExpressionTerm):
    def __init__(self, kernel: Kernel):
        super().__init__(kernel)

    def calculate_coefficients(
        self, point: np.array, boundary_conditions: list[BoundaryCondition], matrix_row: np.array
    ):
        pass

    def calculate_right_side(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        pass

    def value(self, a: float, b: float):
        return a * b


class SingleLayerVolumeTerm(ExpressionTerm):
    def __init__(self, kernel: Kernel):
        super().__init__(kernel)

    def calculate_coefficients(
        self, point: np.array, boundary_conditions: list[BoundaryCondition], matrix_row: np.array
    ):
        pass

    def calculate_right_side(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        pass

    def value(self, a: float, b: float):
        return a * b


print("Define boundary conditions...")

boundary_conditions = [
    BottomDirichletCondition(),
    LeftDirichletCondition(),
    RightDirichletCondition(),
    TopDirichletCondition(),
]

# TODO: Debug
# expression = [SingleLayerBoundaryTerm(), DoubleLayerBoundaryTerm(), SingleLayerVolumeTerm()]

expression = [SingleLayerBoundaryTerm()]

unknown_count = 10

print("Create SLAE...")

matrix = None
right_side = None

for boundary in boundary_conditions:
    for point in boundary.get_boundary_points():
        matrix_row = np.zeros(unknown_count)
        right_side_value = 0.0
        for term in expression:
            term.calculate_coefficients(boundary_conditions, matrix_row)
            right_side_value += term.calculate_right_side()
        if matrix is None:
            matrix = matrix_row
        else:
            matrix = np.vstack((matrix, matrix_row))
        if right_side is None:
            right_side = np.array([right_side_value])
        else:
            right_side = np.append(right_side, [right_side_value])


print("Solving SLAE...")

solution = np.linalg.solve(matrix, right_side)

print("Setting data back...")

for term in expression:
    term.propagate_solution(solution)

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
