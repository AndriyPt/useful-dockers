#!/usr/bin/env python3

from enum import Enum
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

CHART_STEPS = 25
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
        self.value = 0.0


class Utils:
    @staticmethod
    def distance(point1: np.array, point2: np.array):
        return np.linalg.norm(point1 - point2)

    @staticmethod
    def squared_distance(point1: np.array, point2: np.array):
        temp = point1 - point2
        return np.dot(temp.T, temp)

    @staticmethod
    def normalize(vector: np.array):
        norm = np.linalg.norm(vector)
        if 0 == norm:
            vector /= norm
        return vector

    @staticmethod
    def cross_product(vector1: np.array, vector2: np.array):
        if 2 == len(vector1) and 2 == len(vector2):
            value = vector1[0] * vector2[1] - vector1[1] * vector2[0]
            result = np.array([0.0, 0.0, value])
        else:
            result = np.cross(vector1, vector2)
        return result

    @staticmethod
    def is_point_within_segment(point: np.array, begin: np.array, end: np.array, eps: float):
        cross_product = Utils.cross_product(point - begin, end - begin)
        return np.linalg.norm(cross_product) < eps


class BoundaryCondition:
    def __init__(self, border_elements: np.array, type: BoundaryConditionType):
        assert border_elements is not None
        self.__points = []
        for index in range(len(border_elements) - 1):
            item = Point2DInfo()
            item.point = (border_elements[index + 1] + border_elements[index]) / 2.0
            item.type = type

            border_vector = border_elements[index + 1] - border_elements[index]
            assert np.linalg.norm(border_vector) != 0, "Points of border element should be different"
            border_vector /= np.linalg.norm(border_vector)

            # TODO: Check if it is always external normal
            item.normal = np.array([border_vector[1], -border_vector[0]])
            item.element = np.array([border_elements[index], border_elements[index + 1]])
            item.value = self.value(item.point)

            self.__points.append(item)

    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def get_border_elements(self):
        return np.array([item.element[0] for item in self.__points] + [self.__points[-1].element[1]])

    def get_boundary_points(self):
        return [item.point for item in self.__points]

    def get_boundary_points_info(self):
        return self.__points


class DirichletBoundaryCondition(BoundaryCondition):
    def __init__(self, border_elements: np.array):
        super().__init__(border_elements, BoundaryConditionType.DIRICHLET)

    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")


class BottomDirichletCondition(DirichletBoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([0.0, 0.0]), np.array([1.0, 0.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 0.0


class LeftDirichletCondition(DirichletBoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([0.0, 0.0]), np.array([0.0, 1.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 2.0 * point[1]


class RightDirichletCondition(DirichletBoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([1.0, 0.0]), np.array([1.0, 1.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 2.0 * point[1]


class TopDirichletCondition(DirichletBoundaryCondition):
    def __init__(self):
        super().__init__(np.linspace(np.array([0.0, 1.0]), np.array([1.0, 1.0]), BORDER_ELEMENTS_COUNT))

    def value(self, point: np.array):
        return 2.0


class Kernel:
    def value(self, point_x: np.array, point_y: np.array, normal: np.array):
        raise NotImplementedError("Call to abstract method")


class Laplace2DKernel:
    def value(self, point_x: np.array, point_y: np.array, normal: np.array):
        return 0.25 / np.pi * np.log(Utils.squared_distance(point_x, point_y))


class Laplace2DNormKernel:
    def value(self, point_x: np.array, point_y: np.array, normal: np.array):
        raise NotImplementedError("Implement")


class Integrator:
    def __init__(self):
        pass

    def calculate(self, kernel: Kernel, normal: np.array, point: np.array, min_limit: np.array, max_limit: np.array):
        raise NotImplementedError("Call to abstract method")

    def convert_leggauss_to_segment(self, nodes: np.array, begin: np.array, end: np.array):
        result = [(node - (-1.0)) * (end - begin) / 2.0 + begin for node in nodes]
        return result


class Integrator1D(Integrator):

    epsilon = 0.01

    def __init__(self, nominal_count: int, singularity_count: int):
        super().__init__()

        assert 0 < nominal_count and nominal_count < 10, "Not supported points count"
        assert nominal_count < singularity_count and singularity_count < 10, "Not supported points count"

        self.__nominal_number_of_points = nominal_count
        self.__singular_number_of_points = singularity_count

    def calculate(self, kernel: Kernel, normal: np.array, point: np.array, min_limit: np.array, max_limit: np.array):
        result = 0.0
        count = self.__nominal_number_of_points
        if Utils.is_point_within_segment(point, min_limit, max_limit, Integrator1D.epsilon):
            count = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(count)
        real_nodes = self.convert_leggauss_to_segment(nodes, min_limit, max_limit)

        for real_node, weight in zip(real_nodes, weights):
            result += weight * kernel.value(point, real_node, normal)
        return result


class ExpressionTerm:
    def __init__(self, kernel: Kernel):
        self.__kernel = kernel
        assert self.__kernel is not None

    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def kernel(self):
        return self.__kernel

    def calculate_coefficients(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        raise NotImplementedError("Call to abstract method")

    def propagate_solution(self, solution: np.array):
        raise NotImplementedError("Call to abstract method")


class SingleLayerBoundaryTerm(ExpressionTerm):

    NOMINAL_INTEGRATION_POINTS = 4
    SINGULARITY_INTEGRATION_POINTS = 6

    def __init__(self, kernel: Kernel):
        super().__init__(kernel)
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)

    def calculate_coefficients(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        integrator = Integrator1D(
            SingleLayerBoundaryTerm.NOMINAL_INTEGRATION_POINTS, SingleLayerBoundaryTerm.SINGULARITY_INTEGRATION_POINTS
        )
        coefficients = np.empty(0)
        right_side_ret = 0.0

        for boundary_condition in boundary_conditions:
            for boundary_item in boundary_condition.get_boundary_points_info():
                if BoundaryConditionType.DIRICHLET == boundary_item.type:
                    res = integrator.calculate(
                        self.kernel(), boundary_item.normal, point, boundary_item.element[0], boundary_item.element[1]
                    )
                    coefficients = np.append(coefficients, [res])
                elif BoundaryConditionType.NEUMANN == boundary_item.type:
                    right_side_ret += integrator.calculate(
                        self.kernel(), boundary_item.normal, point, boundary_item.point
                    )
                elif BoundaryConditionType.ROBIN == boundary_item.type:
                    raise NotImplementedError("Not implemented")
                else:
                    raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (coefficients, right_side_ret)

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        assert self.__unknown_count == len(self.__unknown_values), "Data should be calculated"
        integrator = Integrator1D(
            SingleLayerBoundaryTerm.NOMINAL_INTEGRATION_POINTS, SingleLayerBoundaryTerm.SINGULARITY_INTEGRATION_POINTS
        )
        result = 0.0
        unknown_index = 0

        for boundary_condition in boundary_conditions:
            for boundary_item in boundary_condition.get_boundary_points_info():
                integral_value = integrator.calculate(
                    self.kernel(), boundary_item.normal, point, boundary_item.element[0], boundary_item.element[1]
                )
                if BoundaryConditionType.DIRICHLET == boundary_item.type:
                    result += self.__unknown_values[unknown_index] * integral_value
                    unknown_index += 1
                elif BoundaryConditionType.NEUMANN == boundary_item.type:
                    result += boundary_item.value * integral_value
                elif BoundaryConditionType.ROBIN == boundary_item.type:
                    raise NotImplementedError("Not implemented")
                else:
                    raise AttributeError("Not supported boundary element type")

        assert self.__unknown_count == unknown_index, "Unknown could should match {} and {}".format(
            self.__unknown_count, unknown_index
        )
        return result


class DoubleLayerBoundaryTerm(ExpressionTerm):
    def __init__(self, kernel: Kernel):
        super().__init__(kernel)

    def calculate_coefficients(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
        pass

    def value(self, a: float, b: float):
        return a * b


class SingleLayerVolumeTerm(ExpressionTerm):
    def __init__(self, kernel: Kernel):
        super().__init__(kernel)

    def calculate_coefficients(self, point: np.array, boundary_conditions: list[BoundaryCondition]):
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

expression = [SingleLayerBoundaryTerm(Laplace2DKernel())]

unknown_count = 10

print("Create SLAE...")

matrix = None
right_side = None

for boundary in boundary_conditions:
    for point in boundary.get_boundary_points():
        right_side_value = 0.0
        matrix_row = np.empty(0)
        for term in expression:
            coefficients, value = term.calculate_coefficients(point, boundary_conditions)
            right_side_value += value
            matrix_row = np.hstack((matrix_row, coefficients))
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
    solution = term.propagate_solution(solution)

print("Visualizing data...")

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

x_min = min(np.min(boundary_condition.get_border_elements()[:, 0]) for boundary_condition in boundary_conditions)
y_min = min(np.min(boundary_condition.get_border_elements()[:, 1]) for boundary_condition in boundary_conditions)
x_max = max(np.max(boundary_condition.get_border_elements()[:, 0]) for boundary_condition in boundary_conditions)
y_max = max(np.max(boundary_condition.get_border_elements()[:, 1]) for boundary_condition in boundary_conditions)

x_data_linear = np.linspace(x_min, x_max, CHART_STEPS)
y_data_linear = np.linspace(y_min, y_max, CHART_STEPS)

x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)

z_data = np.empty([CHART_STEPS, CHART_STEPS])
for x_index in range(CHART_STEPS):
    for y_index in range(CHART_STEPS):
        point = np.array([x_data_linear[x_index], y_data_linear[y_index]])
        z_data[x_index][y_index] = sum(term.value(point, boundary_conditions) for term in expression)

# TODO: Work on numpy way of data visualization
# z_data = sum(np.vectorize(term.value)(x_data, y_data) for term in expression)

surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0)

plt.show()


print("Done!")
