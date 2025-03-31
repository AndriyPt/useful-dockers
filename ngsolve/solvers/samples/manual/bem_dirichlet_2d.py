#!/usr/bin/env python3

from enum import Enum
import itertools
from collections.abc import Callable
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

    @staticmethod
    def is_point_within_square(point: np.array, mesh: np.array, eps: float):
        min_x = mesh[0][0] + eps
        max_x = mesh[1][0] - eps
        min_y = mesh[0][1] + eps
        max_y = mesh[2][1] - eps
        if point[0] > min_x and point[0] < max_x and point[1] > min_y and point[1] < max_y:
            return True
        return False


class Domain:
    def __init__(self):
        pass

    def get_border(self):
        raise NotImplementedError("Call to abstract method")

    def get_mesh(self):
        raise NotImplementedError("Call to abstract method")


class Domain2D(Domain):
    def __init__(self):
        pass

    @staticmethod
    def _process_points(border_elements: np.array, type: BoundaryConditionType, value_function: Callable):
        assert border_elements is not None
        assert value_function is not None
        result = []
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
            item.value = value_function(item.point)

            result.append(item)
        return result

    def get_border(self):
        raise NotImplementedError("Call to abstract method")

    def get_mesh(self):
        raise NotImplementedError("Call to abstract method")


class SquareDomain2D(Domain2D):
    def __init__(
        self,
        bottom_left_point: np.array,
        top_right_point: np.array,
        conditions: list[BoundaryConditionType],
        values: list[Callable],
        side_elements_count: int,
    ):
        super().__init__()

        assert bottom_left_point is not None
        assert top_right_point is not None
        assert 4 == len(conditions)
        assert 4 == len(values)
        assert side_elements_count > 0

        bottom_right_point = np.array([top_right_point[0], bottom_left_point[1]])
        top_left_point = np.array([bottom_left_point[0], top_right_point[1]])

        bottom_points = np.linspace(bottom_left_point, bottom_right_point, side_elements_count)
        left_points = np.linspace(top_left_point, bottom_left_point, side_elements_count)

        self.__border = np.array(Domain2D._process_points(bottom_points, conditions[0], values[0]))
        self.__border = np.append(
            self.__border,
            Domain2D._process_points(
                np.linspace(bottom_right_point, top_right_point, side_elements_count), conditions[1], values[1]
            ),
        )
        self.__border = np.append(
            self.__border,
            Domain2D._process_points(
                np.linspace(top_right_point, top_left_point, side_elements_count), conditions[2], values[2]
            ),
        )
        self.__border = np.array(Domain2D._process_points(left_points, conditions[3], values[3]))

        x_points, y_points = np.meshgrid(bottom_points, np.flip(left_points))

        self.__mesh = []
        for x_index in range(side_elements_count - 1):
            for y_index in range(side_elements_count - 1):
                square = [
                    np.array(x_points[x_index][y_index], y_points[x_index][y_index]),
                    np.array(x_points[x_index + 1][y_index], y_points[x_index + 1][y_index]),
                    np.array(x_points[x_index + 1][y_index + 1], y_points[x_index + 1][y_index + 1]),
                    np.array(x_points[x_index][y_index + 1], y_points[x_index][y_index + 1]),
                ]
                self.__mesh.append(square)
        self.__mesh = np.array(self.__mesh)

    def get_border(self):
        return self.__border

    def get_mesh(self):
        return self.__mesh


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

    def segment(self, kernel: Kernel, normal: np.array, point: np.array, min_limit: np.array, max_limit: np.array):
        raise NotImplementedError("Call to abstract method")

    def square(self, kernel: Kernel, function: Callable, point: np.array, square: np.array):
        raise NotImplementedError("Call to abstract method")


class Integrator2D(Integrator):

    epsilon = 0.01

    def __init__(self, nominal_count: int, singularity_count: int):
        super().__init__()

        assert 0 < nominal_count and nominal_count < 10, "Not supported points count"
        assert nominal_count < singularity_count and singularity_count < 10, "Not supported points count"

        self.__nominal_number_of_points = nominal_count
        self.__singular_number_of_points = singularity_count

    def _convert_leggauss_to_segment(self, nodes: np.array, begin: np.array, end: np.array):
        result = [(node - (-1.0)) * (end - begin) / 2.0 + begin for node in nodes]
        return result

    def _convert_leggauss_to_square(self, nodes: np.array, weights: np.array, square: np.array):

        # TODO: Implement
        return result

    def segment(self, kernel: Kernel, normal: np.array, point: np.array, min_limit: np.array, max_limit: np.array):
        result = 0.0
        count = self.__nominal_number_of_points
        if Utils.is_point_within_segment(point, min_limit, max_limit, Integrator2D.epsilon):
            count = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(count)
        real_nodes = self._convert_leggauss_to_segment(nodes, min_limit, max_limit)

        for real_node, weight in zip(real_nodes, weights):
            result += weight * kernel.value(point, real_node, normal)
        return result

    def square(self, kernel: Kernel, function: Callable, point: np.array, square: np.array):
        result = 0.0
        count = self.__nominal_number_of_points
        if Utils.is_point_within_square(point, square, Integrator2D.epsilon):
            count = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(count)
        (real_nodes, real_weights) = self._convert_leggauss_to_square(nodes, weights, square)

        for real_node, weight in zip(real_nodes, real_weights):
            result += weight * kernel.value(point, real_node, np.empty(2)) * function(real_node)
        return result


class ExpressionTerm:
    def __init__(self, kernel: Kernel):
        self.__kernel = kernel
        assert self.__kernel is not None

    def value(self, point: np.array, domain: Domain):
        raise NotImplementedError("Call to abstract method")

    def kernel(self):
        return self.__kernel

    def calculate_coefficients(self, point: np.array, domain: Domain):
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

    def calculate_coefficients(self, point: np.array, domain: Domain):
        integrator = Integrator2D(
            SingleLayerBoundaryTerm.NOMINAL_INTEGRATION_POINTS, SingleLayerBoundaryTerm.SINGULARITY_INTEGRATION_POINTS
        )
        coefficients = np.empty(0)
        right_side_ret = 0.0

        for boundary_item in domain.get_border():
            if BoundaryConditionType.DIRICHLET == boundary_item.type:
                res = integrator.segment(
                    self.kernel(), boundary_item.normal, point, boundary_item.element[0], boundary_item.element[1]
                )
                coefficients = np.append(coefficients, [res])
            elif BoundaryConditionType.NEUMANN == boundary_item.type:
                right_side_ret += integrator.segment(self.kernel(), boundary_item.normal, point, boundary_item.point)
            elif BoundaryConditionType.ROBIN == boundary_item.type:
                raise NotImplementedError("Not implemented")
            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (coefficients, right_side_ret)

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array, domain: Domain):
        assert self.__unknown_count == len(self.__unknown_values), "Data should be calculated"
        integrator = Integrator2D(
            SingleLayerBoundaryTerm.NOMINAL_INTEGRATION_POINTS, SingleLayerBoundaryTerm.SINGULARITY_INTEGRATION_POINTS
        )
        result = 0.0
        unknown_index = 0

        for boundary_item in domain.get_border():
            integral_value = integrator.segment(
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

    def calculate_coefficients(self, point: np.array, domain: Domain):
        pass

    def value(self, point: np.array, domain: Domain):
        return 0.0


class SingleLayerVolumeTerm(ExpressionTerm):

    NOMINAL_INTEGRATION_POINTS_PER_AXIS = 4
    SINGULARITY_INTEGRATION_POINTS_PER_AXIS = 6

    def __init__(self, kernel: Kernel, value_function: Callable):
        super().__init__(kernel)
        assert value_function is not None
        self.__value_function = value_function

    def calculate_coefficients(self, point: np.array, domain: Domain):
        assert point is not None
        assert domain is not None

        coefficients = np.empty(0)
        right_side_ret = 0.0

        for mesh_item in domain.get_mesh():
            right_side_ret += self.__integrator.calculate_square(self.kernel, self.__value_function, point, mesh_item)

        return (coefficients, right_side_ret)

    def propagate_solution(self, solution: np.array):
        return solution

    def value(self, point: np.array, domain: Domain):
        assert point is not None
        assert domain is not None

        integrator = Integrator2D(
            SingleLayerVolumeTerm.NOMINAL_INTEGRATION_POINTS, SingleLayerVolumeTerm.SINGULARITY_INTEGRATION_POINTS
        )
        result = 0.0
        for mesh_item in domain.get_mesh():
            result += integrator.square(self.kernel, self.__value_function, point, mesh_item)
        return result


print("Define problem...")

print("Define boundary conditions...")


def boundary_value(point: np.array):
    return 2.0 * point[1]


domain = SquareDomain2D(
    np.array([0.0, 0.0]),
    np.array([1.0, 1.0]),
    [BoundaryConditionType.DIRICHLET] * 4,
    [boundary_value] * 4,
    BORDER_ELEMENTS_COUNT,
)

print("Define right hand side value...")


def right_hand_side_function(point: np.array):
    return 2.0


# TODO: Debug
# expression = [SingleLayerBoundaryTerm(), DoubleLayerBoundaryTerm(), SingleLayerVolumeTerm()]

expression = [
    SingleLayerBoundaryTerm(Laplace2DKernel()),
    SingleLayerVolumeTerm(Laplace2DKernel(), right_hand_side_function),
]

print("Create SLAE...")

matrix = None
right_side = None

for point_info in domain.get_border():
    right_side_value = 0.0
    matrix_row = np.empty(0)
    for term in expression:
        coefficients, value = term.calculate_coefficients(point, domain)
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

x_min = min(point_info.point[0] for point_info in domain.get_border())
y_min = min(point_info.point[1] for point_info in domain.get_border())
x_max = max(point_info.point[0] for point_info in domain.get_border())
y_max = max(point_info.point[1] for point_info in domain.get_border())

x_data_linear = np.linspace(x_min, x_max, CHART_STEPS)
y_data_linear = np.linspace(y_min, y_max, CHART_STEPS)

x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)

z_data = np.empty([CHART_STEPS, CHART_STEPS])
for x_index in range(CHART_STEPS):
    for y_index in range(CHART_STEPS):
        point = np.array([x_data_linear[x_index], y_data_linear[y_index]])
        z_data[x_index][y_index] = sum(term.value(point, domain) for term in expression)

# TODO: Work on numpy way of data visualization
# z_data = sum(np.vectorize(term.value)(x_data, y_data) for term in expression)

surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0)

plt.show()


print("Done!")
