#!/usr/bin/env python3

from enum import Enum
from collections.abc import Callable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


class GlobalSettings(object):
    CHART_STEPS = 25
    BORDER_ELEMENTS_COUNT = 10
    PLOT_ERROR = False
    COBORDER_DEPTH = 1.0
    """
        1 - Dirichlet BEM, 2 - Neumann BEM, 3 - Robin BEM, 4 - Single Inclusion Dirichlet BEM
        5 - Dirichlet CoBEM, 6 - Neumann CoBEM, 7 - Robin CoBEM, 8 - Single Inclusion Dirichlet CoBEM    
    """
    EXAMPLE_TYPE = 6


class ExpressionTerm:
    pass


class BoundaryConditionType(Enum):
    UNKNOWN = 1
    DIRICHLET = 2
    NEUMANN = 3
    ROBIN = 4
    INCLUSION = 5


class ProblemSolverType(Enum):
    BEM = 1
    COBEM = 2


class Point2DInfo:
    def __init__(self):
        self.point = np.zeros(2)
        self.type = BoundaryConditionType.UNKNOWN
        self.normal = np.zeros(2)
        self.element = np.empty(0)
        self.value = 0.0
        self.robin_coeff = 0.0  # Robin condition is represented as 1 * q = robin_coeff * u + value


class Utils:
    @staticmethod
    def distance(point1: np.array, point2: np.array):
        return np.linalg.norm(point1 - point2)

    @staticmethod
    def squared_distance(point1: np.array, point2: np.array):
        temp = point1 - point2
        return np.dot(temp.T, temp)

    @staticmethod
    def is_the_same_point(point1: np.array, point2: np.array, eps: float):
        result = Utils.distance(point1, point2) < eps
        return result

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
        result = False
        if Utils.is_the_same_point(point, begin, eps) or Utils.is_the_same_point(point, end, eps):
            result = True
        else:
            cross_product = Utils.cross_product(point - begin, end - begin)
            if np.linalg.norm(cross_product) < eps:
                result = True
        return result

    @staticmethod
    def is_point_within_square(point: np.array, mesh: np.array, eps: float):
        min_x = mesh[0][0] + eps
        max_x = mesh[1][0] - eps
        min_y = mesh[0][1] + eps
        max_y = mesh[2][1] - eps
        if point[0] > min_x and point[0] < max_x and point[1] > min_y and point[1] < max_y:
            return True
        return False

    # TODO: Implement more precise logic
    @staticmethod
    def is_point_within_trapezoid(point: np.array, mesh: np.array, eps: float):
        return Utils.is_point_within_square(point, mesh, eps)

    @staticmethod
    def is_square_side(begin: np.array, end: np.array, mesh_item: np.array, eps: float):
        assert 2 == len(begin)
        assert 2 == len(end)
        assert 4 == len(mesh_item)

        result = False
        for index in range(0, 2):
            if (
                Utils.is_the_same_point(begin, mesh_item[index], eps)
                and Utils.is_the_same_point(end, mesh_item[index + 1], eps)
            ) or (
                Utils.is_the_same_point(end, mesh_item[index], eps)
                and Utils.is_the_same_point(begin, mesh_item[index + 1], eps)
            ):
                result = True
                break
        return result

    # TODO: Fix a name of a geometrical shape :)
    @staticmethod
    def is_rombus(points: np.array, eps: float):
        assert 4 == len(points)

        distance = Utils.distance(points[0], points[1])
        distance += Utils.distance(points[1], points[2])
        distance -= Utils.distance(points[2], points[3])
        distance -= Utils.distance(points[3], points[0])

        result = np.abs(distance) < eps
        return result

    @staticmethod
    def constant_one():
        def result(point: np.array):
            return 1.0

        return result

    @staticmethod
    def constant_value(value: np.array):
        def result(point: np.array):
            return value

        return result


class Domain:
    def __init__(self, subdomains: list = []):
        self.__subdomains = subdomains

    def get_border(self):
        raise NotImplementedError("Call to abstract method")

    def get_mesh(self):
        raise NotImplementedError("Call to abstract method")

    def get_coborder(self):
        raise NotImplementedError("Call to abstract method")

    def is_point_on_border(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def get_subdomains(self):
        return self.__subdomains


class Domain2D(Domain):
    POINT_LOCATION_EPSILON = 0.001

    def __init__(self, subdomains: list = []):
        super().__init__(subdomains)

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
            if BoundaryConditionType.ROBIN == type:
                item.robin_coeff, item.value = value_function(item.point)
            else:
                item.value = value_function(item.point)

            result.append(item)
        return result

    def get_border(self):
        raise NotImplementedError("Call to abstract method")

    def get_mesh(self):
        raise NotImplementedError("Call to abstract method")

    def get_coborder(self):
        raise NotImplementedError("Call to abstract method")

    def is_point_on_border(self, point: np.array):
        for border_point in self.get_border():
            if Utils.is_point_within_segment(
                point, border_point.element[0], border_point.element[1], Domain2D.POINT_LOCATION_EPSILON
            ):
                return True
        return False


class SquareDomain2D(Domain2D):
    def __init__(
        self,
        bottom_left_point: np.array,
        top_right_point: np.array,
        conditions: list[BoundaryConditionType],
        values: list[Callable],
        side_elements_count: int,
        subdomains: list = [],
    ):
        super().__init__(subdomains)

        assert bottom_left_point is not None
        assert top_right_point is not None
        assert 4 == len(conditions)
        assert 4 == len(values)
        assert side_elements_count > 0

        side_points_count = side_elements_count + 1

        bottom_right_point = np.array([top_right_point[0], bottom_left_point[1]])
        top_left_point = np.array([bottom_left_point[0], top_right_point[1]])

        bottom_points = np.linspace(bottom_left_point, bottom_right_point, side_points_count)
        left_points = np.linspace(top_left_point, bottom_left_point, side_points_count)

        self.__border = np.array(Domain2D._process_points(bottom_points, conditions[0], values[0]))
        self.__border = np.append(
            self.__border,
            Domain2D._process_points(
                np.linspace(bottom_right_point, top_right_point, side_points_count), conditions[1], values[1]
            ),
        )
        self.__border = np.append(
            self.__border,
            Domain2D._process_points(
                np.linspace(top_right_point, top_left_point, side_points_count), conditions[2], values[2]
            ),
        )
        self.__border = np.append(
            self.__border, np.array(Domain2D._process_points(left_points, conditions[3], values[3]))
        )

        x_points, y_points = np.meshgrid(bottom_points[:, 0], np.flip(left_points[:, 1]))

        self.__mesh = []
        for x_index in range(side_points_count - 1):
            for y_index in range(side_points_count - 1):
                square = [
                    np.array([x_points[x_index][y_index], y_points[x_index][y_index]]),
                    np.array([x_points[x_index][y_index + 1], y_points[x_index][y_index + 1]]),
                    np.array([x_points[x_index + 1][y_index + 1], y_points[x_index + 1][y_index + 1]]),
                    np.array([x_points[x_index + 1][y_index], y_points[x_index + 1][y_index]]),
                ]
                point_info = Point2DInfo()
                point_info.point = (
                    np.array(
                        [
                            x_points[x_index][y_index] + x_points[x_index + 1][y_index + 1],
                            y_points[x_index][y_index] + y_points[x_index + 1][y_index + 1],
                        ]
                    )
                    / 2.0
                )
                point_info.type = BoundaryConditionType.INCLUSION
                point_info.value = 0.0
                point_info.element = np.array(square)
                self.__mesh.append(point_info)
        self.__mesh = np.array(self.__mesh)
        self.__square = np.array([bottom_left_point, bottom_right_point, top_right_point, top_left_point])
        self._fill_coborder_elements()

    def _fill_coborder_elements(self):
        assert self.__border is not None
        assert self.__square is not None

        self.__coborder = []
        for index in range(len(self.__border)):
            point_info = self.__border[index]
            item = Point2DInfo()
            item.point = point_info.point
            item.type = point_info.type
            item.normal = point_info.normal
            item.value = point_info.value
            item.robin_coeff = point_info.robin_coeff

            corner_points = []
            for corner in self.__square:
                for element in point_info.element:
                    if Utils.is_the_same_point(element, corner, Domain2D.POINT_LOCATION_EPSILON):
                        corner_points.append(element)

            elements = []
            elements.append(point_info.element[0])
            elements.append(point_info.element[1])

            if len(corner_points) > 0:
                corner_point = corner_points[0]
                if Utils.is_the_same_point(elements[0], corner_point, Domain2D.POINT_LOCATION_EPSILON) or 2 == len(
                    corner_points
                ):
                    prev_index = index - 1
                    if prev_index < 0:
                        prev_index = len(self.__border) - 1
                    corner_grow_vector = self.__border[prev_index].normal + point_info.normal
                    elements.append(point_info.element[0] + corner_grow_vector * GlobalSettings.COBORDER_DEPTH)
                    if 1 == len(corner_points):
                        elements.append(point_info.element[1] + point_info.normal * GlobalSettings.COBORDER_DEPTH)

                if Utils.is_the_same_point(elements[1], corner_point, Domain2D.POINT_LOCATION_EPSILON) or 2 == len(
                    corner_points
                ):
                    next_index = index + 1
                    if next_index >= len(self.__border):
                        next_index = 0
                    corner_grow_vector = self.__border[next_index].normal + point_info.normal
                    elements.append(point_info.element[1] + corner_grow_vector * GlobalSettings.COBORDER_DEPTH)
                    if 1 == len(corner_points):
                        elements.append(point_info.element[0] + point_info.normal * GlobalSettings.COBORDER_DEPTH)
            else:
                elements.append(point_info.element[1] + point_info.normal * GlobalSettings.COBORDER_DEPTH)
                elements.append(point_info.element[0] + point_info.normal * GlobalSettings.COBORDER_DEPTH)

            if Utils.is_the_same_point(point_info.normal, np.array([1.0, 0.0]), Domain2D.POINT_LOCATION_EPSILON):
                elements[1], elements[2], elements[3] = elements[2], elements[3], elements[1]
            elif Utils.is_the_same_point(point_info.normal, np.array([0.0, 1.0]), Domain2D.POINT_LOCATION_EPSILON):
                elements[0], elements[1] = elements[1], elements[0]
            elif Utils.is_the_same_point(point_info.normal, np.array([-1.0, 0.0]), Domain2D.POINT_LOCATION_EPSILON):
                elements[0], elements[2], elements[3] = elements[3], elements[0], elements[2]
            else:
                elements[3], elements[2], elements[0], elements[1] = elements[0], elements[1], elements[2], elements[3]

            item.element = elements
            self.__coborder.append(item)

        self.__coborder = np.array(self.__coborder)

    def get_border(self):
        return self.__border

    def get_mesh(self):
        return self.__mesh

    def get_coborder(self):
        return self.__coborder

    def get_square(self):
        return self.__square

    def is_point_inside_domain(self, point: np.array):
        result = Utils.is_point_within_square(point, self.__square, Domain2D.POINT_LOCATION_EPSILON)
        return result


class Kernel:
    def value(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def grad(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")


class Laplace2DKernel(Kernel):
    def value(self, point_x: np.array, point_y: np.array):
        return -0.25 / np.pi * np.log(Utils.squared_distance(point_x, point_y))

    def grad(self, point_x: np.array, point_y: np.array):
        result = -0.25 / np.pi / Utils.squared_distance(point_x, point_y) * -2.0 * (point_x - point_y)
        return result


class Integrator:
    def __init__(self, kernel: Kernel):
        assert kernel is not None
        self.__kernel = kernel

    @property
    def kernel(self):
        return self.__kernel

    def segment(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        raise NotImplementedError("Call to abstract method")

    def segment_grad(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        raise NotImplementedError("Call to abstract method")

    def square(self, function: Callable, point: np.array, square: np.array):
        raise NotImplementedError("Call to abstract method")

    def square_grad(self, function: Callable, point: np.array, square: np.array):
        raise NotImplementedError("Call to abstract method")


class Integrator2D(Integrator):

    EPSILON = 0.01
    DIM = 2

    def __init__(self, kernel: Kernel, nominal_count: int, singularity_count: int):
        super().__init__(kernel)

        assert 0 < nominal_count and nominal_count < 10, "Not supported points count"
        assert nominal_count < singularity_count and singularity_count < 10, "Not supported points count"

        self.__nominal_number_of_points = nominal_count
        self.__singular_number_of_points = singularity_count

    def _convert_leggauss_to_segment(self, nodes: np.array, weights: np.array, begin: np.array, end: np.array):
        result_nodes = np.array([(node - (-1.0)) * (end - begin) / 2.0 + begin for node in nodes])
        weights_adjustment = Utils.distance(begin, end) / 2.0
        result_weights = weights_adjustment * weights
        return (result_nodes, result_weights)

    def _convert_leggauss_to_square(self, nodes: np.array, weights: np.array, square: np.array):
        assert 4 == len(square)
        x_nodes, x_weights = self._convert_leggauss_to_segment(nodes, weights, square[0], square[1])
        x_nodes = x_nodes[:, 0]
        y_nodes, y_weights = self._convert_leggauss_to_segment(nodes, weights, square[0], square[3])
        y_nodes = y_nodes[:, 1]

        x_nodes, y_nodes = np.meshgrid(x_nodes, y_nodes)
        x_nodes = np.ndarray.flatten(x_nodes)
        y_nodes = np.ndarray.flatten(y_nodes)
        result_nodes = np.column_stack((x_nodes, y_nodes))

        x_weights, y_weights = np.meshgrid(x_weights, y_weights)
        x_weights = np.ndarray.flatten(x_weights)
        y_weights = np.ndarray.flatten(y_weights)
        result_weights = x_weights * y_weights

        return (result_nodes, result_weights)

    def _get_segment_nodes_and_weights(self, point: np.array, min_limit: np.array, max_limit: np.array):
        count = self.__nominal_number_of_points
        # TODO: Implement logic which will increase points count in case if one point is too close
        if Utils.is_point_within_segment(point, min_limit, max_limit, Integrator2D.EPSILON):
            count = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(count)
        real_nodes, real_weights = self._convert_leggauss_to_segment(nodes, weights, min_limit, max_limit)
        return (real_nodes, real_weights)

    def segment(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        result = 0.0
        real_nodes, real_weights = self._get_segment_nodes_and_weights(point, min_limit, max_limit)
        for real_node, weight in zip(real_nodes, real_weights):
            result += weight * self.kernel.value(point, real_node) * function(point)
        return result

    def segment_grad(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        result = 0.0
        real_nodes, real_weights = self._get_segment_nodes_and_weights(point, min_limit, max_limit)
        for real_node, weight in zip(real_nodes, real_weights):
            result += weight * np.dot(self.kernel.grad(point, real_node), function(point))
        return result

    def _get_square_nodes_and_weights(self, point: np.array, square: np.array):
        count = self.__nominal_number_of_points
        if Utils.is_point_within_square(point, square, Integrator2D.EPSILON):
            count = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(count)
        (real_nodes, real_weights) = self._convert_leggauss_to_square(nodes, weights, square)
        return (real_nodes, real_weights)

    def square(self, function: Callable, point: np.array, square: np.array):
        result = 0.0
        real_nodes, real_weights = self._get_square_nodes_and_weights(point, square)
        for real_node, weight in zip(real_nodes, real_weights):
            result += weight * self.kernel.value(point, real_node) * function(real_node)
        return result

    def square_grad(self, function: Callable, point: np.array, square: np.array):
        result = 0.0
        real_nodes, real_weights = self._get_square_nodes_and_weights(point, square)
        for real_node, weight in zip(real_nodes, real_weights):
            result += weight * np.dot(self.kernel.grad(point, real_node), function(real_node))
        return result

    # Use Jacobian to transform trapezoid to unit square
    # Find transformation between coordinate systems using representation x' = a*x + b*y + c*x*y + const_x
    # and y' = d*x + e*y + f*x*y + const_y
    # Jacobian (a*e - b*d) + (c*e - b*f)*y + (a*f - c*d)*y
    def _get_trapezoid_nodes_weights_and_coeff(self, point: np.array, trapezoid: np.array):
        unit_square = np.array(
            [np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([1.0, 1.0]), np.array([0.0, 1.0])]
        )
        matrix = None
        right_side = np.empty(0)
        for ideal, actual in zip(unit_square, trapezoid):
            matrix_row = np.zeros(2 * len(unit_square))
            matrix_row[0] = 1.0
            matrix_row[1] = ideal[0]
            matrix_row[2] = ideal[1]
            matrix_row[3] = ideal[0] * ideal[1]
            right_side = np.append(right_side, [actual[0]])
            if matrix is None:
                matrix = matrix_row
            else:
                matrix = np.vstack((matrix, matrix_row))
            matrix_row = np.zeros(2 * len(unit_square))
            shift = len(unit_square)
            matrix_row[shift + 0] = 1.0
            matrix_row[shift + 1] = ideal[0]
            matrix_row[shift + 2] = ideal[1]
            matrix_row[shift + 3] = ideal[0] * ideal[1]
            right_side = np.append(right_side, [actual[1]])
            matrix = np.vstack((matrix, matrix_row))

        coeff = np.linalg.solve(matrix, right_side)
        if Utils.is_point_within_trapezoid(point, trapezoid, Integrator2D.EPSILON):
            nodes, weights = self._get_square_nodes_and_weights(np.array([0.5, 0.5]), unit_square)
        else:
            nodes, weights = self._get_square_nodes_and_weights(np.array([10.0, 10.0]), unit_square)
        return (nodes, weights, coeff)

    def trapezoid(self, function: Callable, point: np.array, trapezoid: np.array):
        assert 4 == len(trapezoid)
        result = 0.0
        if Utils.is_rombus(trapezoid, Integrator2D.EPSILON):
            result = self.square(function, point, trapezoid)
            return result
        nodes, weights, coeff = self._get_trapezoid_nodes_weights_and_coeff(point, trapezoid)
        for node, weight in zip(nodes, weights):
            real_node = np.zeros(Integrator2D.DIM)
            a, b, c, d, e, f = (coeff[1], coeff[2], coeff[3], coeff[5], coeff[6], coeff[7])
            const_x, const_y = (coeff[0], coeff[4])
            real_node[0] = const_x + a * node[0] + b * node[1] + c * node[0] * node[1]
            real_node[1] = const_y + d * node[0] + e * node[1] + f * node[0] * node[1]
            jacobian = (a * e - b * d) + (a * f - c * d) * node[0] + (c * e - b * f) * node[1]
            result += weight * self.kernel.value(point, real_node) * function(real_node) * jacobian
        return result

    def trapezoid_grad(self, function: Callable, point: np.array, trapezoid: np.array):
        assert 4 == len(trapezoid)
        result = 0.0
        if Utils.is_rombus(trapezoid, Integrator2D.EPSILON):
            result = self.square_grad(function, point, trapezoid)
            return result
        nodes, weights, coeff = self._get_trapezoid_nodes_weights_and_coeff(point, trapezoid)
        for node, weight in zip(nodes, weights):
            real_node = np.zeros(Integrator2D.DIM)
            a, b, c, d, e, f = (coeff[1], coeff[2], coeff[3], coeff[5], coeff[6], coeff[7])
            const_x, const_y = (coeff[0], coeff[4])
            real_node[0] = const_x + a * node[0] + b * node[1] + c * node[0] * node[1]
            real_node[1] = const_y + d * node[0] + e * node[1] + f * node[0] * node[1]
            jacobian = (a * e - b * d) + (a * f - c * d) * node[0] + (c * e - b * f) * node[1]
            result += weight * np.dot(self.kernel.grad(point, real_node), function(real_node)) * jacobian
        return result


class ExpressionTerm:
    def __init__(self, kernel: Kernel, domain: Domain, sign: int):
        assert kernel is not None
        assert domain is not None
        assert 1 == sign or -1 == sign

        self.__kernel = kernel
        self.__domain = domain
        self.__sign = float(sign)

    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    @property
    def kernel(self):
        return self.__kernel

    @property
    def domain(self):
        return self.__domain

    @property
    def sign(self):
        return self.__sign

    def calculate_coefficients(self, point_info: Point2DInfo):
        raise NotImplementedError("Call to abstract method")

    def calculate_for_robin(self, point_info: Point2DInfo):
        raise NotImplementedError("Call to abstract method")

    def propagate_solution(self, solution: np.array):
        raise NotImplementedError("Call to abstract method")


class SingleLayerBoundaryTerm(ExpressionTerm):

    NOMINAL_INTEGRATION_POINTS = 4
    SINGULARITY_INTEGRATION_POINTS = 6
    EPS = 0.001

    def __init__(self, kernel: Kernel, domain: Domain, sign: int = 1):
        super().__init__(kernel, domain, sign)
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)
        self.__integrator = Integrator2D(
            kernel,
            SingleLayerBoundaryTerm.NOMINAL_INTEGRATION_POINTS,
            SingleLayerBoundaryTerm.SINGULARITY_INTEGRATION_POINTS,
        )

    def calculate_coefficients(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0

        for boundary_item in self.domain.get_border():
            if boundary_item.type in [BoundaryConditionType.DIRICHLET, BoundaryConditionType.ROBIN]:
                res = self.__integrator.segment(
                    Utils.constant_one(), point_info.point, boundary_item.element[0], boundary_item.element[1]
                )
                coefficients = np.append(coefficients, [res])
            elif BoundaryConditionType.NEUMANN == boundary_item.type:
                right_side_ret += boundary_item.value * self.__integrator.segment(
                    Utils.constant_one(), point_info.point, boundary_item.element[0], boundary_item.element[1]
                )
            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (self.sign * coefficients, self.sign * right_side_ret)

    def calculate_for_robin(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        for boundary_item in self.domain.get_border():
            mid_point = 0.5 * (boundary_item.element[1] + boundary_item.element[0])
            is_same_point = Utils.distance(point_info.point, mid_point) < SingleLayerBoundaryTerm.EPS

            if BoundaryConditionType.DIRICHLET == boundary_item.type:
                coefficients = np.append(coefficients, [0.0])
            elif BoundaryConditionType.NEUMANN == boundary_item.type:
                pass
            elif BoundaryConditionType.ROBIN == boundary_item.type:
                if is_same_point:
                    coefficients = np.append(coefficients, [1.0])
                else:
                    coefficients = np.append(coefficients, [0.0])
            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (self.sign * coefficients, 0.0)

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array):
        assert self.__unknown_count == len(self.__unknown_values), "Data should be calculated"
        result = 0.0
        unknown_index = 0

        for boundary_item in self.domain.get_border():
            integral_value = self.__integrator.segment(
                Utils.constant_one(), point, boundary_item.element[0], boundary_item.element[1]
            )
            if boundary_item.type in [BoundaryConditionType.DIRICHLET, BoundaryConditionType.ROBIN]:
                result += self.__unknown_values[unknown_index] * integral_value
                unknown_index += 1
            elif BoundaryConditionType.NEUMANN == boundary_item.type:
                result += boundary_item.value * integral_value
            else:
                raise AttributeError("Not supported boundary element type")

        assert self.__unknown_count == unknown_index, "Unknown could should match {} and {}".format(
            self.__unknown_count, unknown_index
        )
        result *= self.sign
        return result


class DoubleLayerBoundaryTerm(ExpressionTerm):
    NOMINAL_INTEGRATION_POINTS = 4
    SINGULARITY_INTEGRATION_POINTS = 6
    EPS = 0.001

    def __init__(self, kernel: Kernel, domain: Domain, sign: int = 1):
        super().__init__(kernel, domain, sign)
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)
        self.__integrator = Integrator2D(
            kernel,
            DoubleLayerBoundaryTerm.NOMINAL_INTEGRATION_POINTS,
            DoubleLayerBoundaryTerm.SINGULARITY_INTEGRATION_POINTS,
        )

    def calculate_coefficients(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0
        for boundary_item in self.domain.get_border():
            is_same_point = Utils.is_the_same_point(point_info.point, boundary_item.point, DoubleLayerBoundaryTerm.EPS)
            if BoundaryConditionType.DIRICHLET == boundary_item.type:
                right_side_ret += boundary_item.value * self.__integrator.segment_grad(
                    Utils.constant_value(boundary_item.normal),
                    point_info.point,
                    boundary_item.element[0],
                    boundary_item.element[1],
                )
                if is_same_point:
                    right_side_ret += 0.5 * boundary_item.value
            elif boundary_item.type in [BoundaryConditionType.NEUMANN, BoundaryConditionType.ROBIN]:
                res = self.__integrator.segment_grad(
                    Utils.constant_value(boundary_item.normal),
                    point_info.point,
                    boundary_item.element[0],
                    boundary_item.element[1],
                )
                if is_same_point:
                    res += 0.5
                coefficients = np.append(coefficients, [res])
            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (self.sign * coefficients, self.sign * right_side_ret)

    def calculate_for_robin(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0
        for boundary_item in self.domain.get_border():
            mid_point = 0.5 * (boundary_item.element[1] + boundary_item.element[0])
            is_same_point = Utils.distance(point_info.point, mid_point) < SingleLayerBoundaryTerm.EPS

            if BoundaryConditionType.DIRICHLET == boundary_item.type:
                pass
            elif BoundaryConditionType.NEUMANN == boundary_item.type:
                coefficients = np.append(coefficients, [0.0])
            elif BoundaryConditionType.ROBIN == boundary_item.type:
                if is_same_point:
                    right_side_ret += boundary_item.value
                    coefficients = np.append(coefficients, [boundary_item.robin_coeff])
                else:
                    coefficients = np.append(coefficients, [0.0])
            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (self.sign * coefficients, self.sign * right_side_ret)

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array):
        assert self.__unknown_count == len(self.__unknown_values), "Data should be calculated"
        result = 0.0
        unknown_index = 0

        for boundary_item in self.domain.get_border():
            integral_value = self.__integrator.segment_grad(
                Utils.constant_value(boundary_item.normal), point, boundary_item.element[0], boundary_item.element[1]
            )
            if BoundaryConditionType.DIRICHLET == boundary_item.type:
                result += boundary_item.value * integral_value
            elif boundary_item.type in [BoundaryConditionType.NEUMANN, BoundaryConditionType.ROBIN]:
                result += self.__unknown_values[unknown_index] * integral_value
                unknown_index += 1
            else:
                raise AttributeError("Not supported boundary element type")

        assert self.__unknown_count == unknown_index, "Unknown could should match {} and {}".format(
            self.__unknown_count, unknown_index
        )
        result *= self.sign
        return result


class SingleLayerVolumeTerm(ExpressionTerm):

    NOMINAL_INTEGRATION_POINTS_PER_AXIS = 2
    SINGULARITY_INTEGRATION_POINTS_PER_AXIS = 4

    def __init__(self, kernel: Kernel, domain: Domain, value_function: Callable, sign: int = 1):
        super().__init__(kernel, domain, sign)
        assert value_function is not None
        self.__value_function = value_function
        self.__integrator = Integrator2D(
            kernel,
            SingleLayerVolumeTerm.NOMINAL_INTEGRATION_POINTS_PER_AXIS,
            SingleLayerVolumeTerm.SINGULARITY_INTEGRATION_POINTS_PER_AXIS,
        )

    def _get_integrator(self):
        return self.__integrator

    def _get_value_function(self):
        return self.__value_function

    def calculate_coefficients(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = self.value(point_info.point)
        return (coefficients, right_side_ret)

    def calculate_for_robin(self, point_info: Point2DInfo):
        coefficients = np.empty(0)
        return (coefficients, 0.0)

    def propagate_solution(self, solution: np.array):
        return solution

    def value(self, point: np.array):
        assert point is not None
        result = 0.0
        for point_info in self.domain.get_mesh():
            result += self.__integrator.square(self.__value_function, point, point_info.element)
        result *= self.sign
        return result


class SingleLayerVolumeCoBEMTerm(SingleLayerVolumeTerm):

    EPS = 0.001

    def __init__(self, kernel: Kernel, domain: Domain, value_function: Callable, sign: int = 1):
        super().__init__(kernel, domain, value_function, sign)

    def calculate_coefficients(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        if BoundaryConditionType.NEUMANN == point_info.type:
            right_side_ret = self.value_grad(point_info.point, point_info.normal)
        else:
            right_side_ret = super().value(point_info.point)
        return (coefficients, right_side_ret)

    def value_grad(self, point: np.array, norm: np.array):
        assert point is not None
        result = 0.0

        def normal_function(x: np.array):
            result = self._get_value_function()(x) * norm
            return result

        for point_info in self.domain.get_mesh():
            result += self._get_integrator().square_grad(normal_function, point, point_info.element)
        result *= self.sign
        return result


class SingleLayerCoBoundaryTerm(ExpressionTerm):
    NOMINAL_INTEGRATION_POINTS = 4
    SINGULARITY_INTEGRATION_POINTS = 6
    EPS = 0.001

    def __init__(self, kernel: Kernel, domain: Domain, sign: int = 1):
        super().__init__(kernel, domain, sign)
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)
        self.__integrator = Integrator2D(
            kernel,
            SingleLayerCoBoundaryTerm.NOMINAL_INTEGRATION_POINTS,
            SingleLayerCoBoundaryTerm.SINGULARITY_INTEGRATION_POINTS,  # There should not be singularities
        )

    def calculate_coefficients(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0
        for boundary_item in self.domain.get_coborder():
            is_same_point = Utils.is_the_same_point(
                point_info.point, boundary_item.point, SingleLayerCoBoundaryTerm.EPS
            )
            if BoundaryConditionType.DIRICHLET == point_info.type:
                res = self.__integrator.trapezoid(Utils.constant_one(), point_info.point, boundary_item.element)
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value
            elif BoundaryConditionType.NEUMANN == point_info.type:
                res = self.__integrator.trapezoid_grad(
                    Utils.constant_value(boundary_item.normal), point_info.point, boundary_item.element
                )
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value
            elif BoundaryConditionType.ROBIN == point_info.type:
                res = self.__integrator.trapezoid(Utils.constant_one(), point_info.point, boundary_item.element)
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value
                res = self.__integrator.trapezoid_grad(
                    Utils.constant_value(boundary_item.normal), point_info.point, boundary_item.element
                )
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value

            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (self.sign * coefficients, self.sign * right_side_ret)

    def calculate_for_robin(self, point_info: Point2DInfo):
        raise NotImplementedError("Not implemented yet")

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array):
        assert self.__unknown_count == len(self.__unknown_values), "Data should be calculated"
        result = 0.0
        unknown_index = 0

        for boundary_item in self.domain.get_coborder():
            integral_value = self.__integrator.trapezoid(Utils.constant_one(), point, boundary_item.element)
            if boundary_item.type in [BoundaryConditionType.DIRICHLET, BoundaryConditionType.NEUMANN]:
                result += self.__unknown_values[unknown_index] * integral_value
                unknown_index += 1
            else:
                raise AttributeError("Not supported boundary element type")

        assert self.__unknown_count == unknown_index, "Unknown could should match {} and {}".format(
            self.__unknown_count, unknown_index
        )
        result *= self.sign
        return result


class SingleLayerInclusionTerm(ExpressionTerm):

    NOMINAL_INTEGRATION_POINTS_PER_AXIS = 2
    SINGULARITY_INTEGRATION_POINTS_PER_AXIS = 4

    def __init__(
        self, kernel: Kernel, domain: Domain, grad_function: Callable, laplacian_function: Callable, sign: int = 1
    ):
        super().__init__(kernel, domain, sign)
        assert grad_function is not None
        assert laplacian_function is not None
        self.__grad_function = grad_function
        self.__laplacian_function = laplacian_function
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)
        self.__integrator = Integrator2D(
            kernel,
            SingleLayerInclusionTerm.NOMINAL_INTEGRATION_POINTS_PER_AXIS,
            SingleLayerInclusionTerm.SINGULARITY_INTEGRATION_POINTS_PER_AXIS,
        )

    def calculate_coefficients(self, point_info: Point2DInfo):
        assert point_info is not None
        coefficients = np.empty(0)
        for point_info in self.domain.get_mesh():
            is_same_point = Utils.is_the_same_point(point_info.point, point_info.point, DoubleLayerBoundaryTerm.EPS)
            res = -1.0 * self.__integrator.square(self.__laplacian_function, point_info.point, point_info.element)
            res += self.__integrator.square_grad(self.__grad_function, point_info.point, point_info.element)
            # TODO: Check sign of the U variable
            if is_same_point:
                res += 1.0
            for boundary_item in self.domain.get_border():
                if Utils.is_square_side(
                    boundary_item.element[0],
                    boundary_item.element[1],
                    point_info.element,
                    Domain2D.POINT_LOCATION_EPSILON,
                ):

                    def normal_derivative(position: np.array):
                        return np.dot(self.__grad_function(position), boundary_item.normal)

                    res += self.__integrator.segment(
                        normal_derivative, point_info.point, boundary_item.element[0], boundary_item.element[1]
                    )

            coefficients = np.append(coefficients, [res])

        self.__unknown_count = len(coefficients)
        return (self.sign * coefficients, 0.0)

    def calculate_for_robin(self, point_info: Point2DInfo):
        raise NotImplementedError("Implement")

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array):
        assert point is not None
        result = 0.0
        unknown_index = 0
        for point_info in self.domain.get_mesh():
            result -= self.__integrator.square(self.__laplacian_function, point, point_info.element)
            result += self.__integrator.square_grad(self.__grad_function, point, point_info.element)
            for boundary_item in self.domain.get_border():
                if Utils.is_square_side(
                    boundary_item.element[0],
                    boundary_item.element[1],
                    point_info.element,
                    Domain2D.POINT_LOCATION_EPSILON,
                ):

                    def normal_derivative(position: np.array):
                        return np.dot(self.__grad_function(position), boundary_item.normal)

                    result += self.__integrator.segment(
                        normal_derivative, point, boundary_item.element[0], boundary_item.element[1]
                    )
            result *= self.__unknown_values[unknown_index]

            unknown_index += 1

        assert self.__unknown_count == unknown_index, "Unknown could should match {} and {}".format(
            self.__unknown_count, unknown_index
        )
        result *= self.sign
        return result


class Problem(object):
    def __init__(
        self,
        type: ProblemSolverType,
        expression: list[ExpressionTerm],
        domain: Domain,
        analytical_solution: Callable = None,
    ):
        assert type is not None
        assert expression is not None
        assert domain is not None
        self.__type = type
        self.__expression = expression
        self.__domain = domain
        self.__analytical_solution = analytical_solution

    def solution_value(self, point: np.array):
        result = -1.0 * sum(term.value(point) for term in self.__expression)
        if ProblemSolverType.BEM == self.__type:
            # TODO: Check this behaviour
            if self.__domain.is_point_on_border(point):
                result *= 2.0
            for corner in self.__domain.get_square():
                if Utils.distance(point, corner) < Domain2D.POINT_LOCATION_EPSILON:
                    result *= 2.0
        return result

    def calculate(self):
        print("Create SLAE...")

        matrix = None
        right_side = None

        point_list = self.__domain.get_border()
        for subdomain in self.__domain.get_subdomains():
            point_list = np.hstack((point_list, subdomain.get_mesh()))

        for point_info in point_list:
            right_side_value = 0.0
            matrix_row = np.empty(0)
            for term in self.__expression:
                coefficients, value = term.calculate_coefficients(point_info)
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
            if BoundaryConditionType.ROBIN == point_info.type and ProblemSolverType.BEM == self.__type:
                right_side_value = 0.0
                matrix_row = np.empty(0)
                for term in self.__expression:
                    coefficients, value = term.calculate_for_robin(point_info)
                    right_side_value += value
                    matrix_row = np.hstack((matrix_row, coefficients))
                matrix = np.vstack((matrix, matrix_row))
                right_side = np.append(right_side, [right_side_value])

        print("Solving SLAE...")

        solution = np.linalg.solve(matrix, -1.0 * right_side)

        print("Setting data back...")

        for term in self.__expression:
            solution = term.propagate_solution(solution)

    def plot(self):
        print("Visualizing data...")

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        x_min = min(point_info.point[0] for point_info in self.__domain.get_border())
        y_min = min(point_info.point[1] for point_info in self.__domain.get_border())
        x_max = max(point_info.point[0] for point_info in self.__domain.get_border())
        y_max = max(point_info.point[1] for point_info in self.__domain.get_border())

        x_data_linear = np.linspace(x_min, x_max, GlobalSettings.CHART_STEPS)
        y_data_linear = np.linspace(y_min, y_max, GlobalSettings.CHART_STEPS)

        x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)

        z_data = np.empty([GlobalSettings.CHART_STEPS, GlobalSettings.CHART_STEPS])
        for x_index in range(GlobalSettings.CHART_STEPS):
            for y_index in range(GlobalSettings.CHART_STEPS):
                point = np.array([x_data_linear[x_index], y_data_linear[y_index]])
                if GlobalSettings.PLOT_ERROR and self.__analytical_solution is not None:
                    z_data[x_index][y_index] = np.abs(self.solution_value(point) - self.__analytical_solution(point))
                else:
                    z_data[x_index][y_index] = self.solution_value(point)

        # TODO: Work on numpy way of data visualization
        # z_data = sum(np.vectorize(term.value)(x_data, y_data) for term in expression)

        surf = ax.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0)

        plt.show()

        print("Done!")


def init_poisson_dirichlet_bem():
    print("BEM for Dirichlet problem for Poisson equation...")

    print("Define boundary conditions...")

    def analytical_solution(point: np.array):
        return 2 * (point[0] - 0.5) ** 2 + 2 * (point[1] - 0.5) ** 2

    def boundary_value(point: np.array):
        return analytical_solution(point)

    domain = SquareDomain2D(
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0]),
        [BoundaryConditionType.DIRICHLET] * 4,
        [boundary_value] * 4,
        GlobalSettings.BORDER_ELEMENTS_COUNT,
    )

    print("Define heat source function...")

    def heat_source_function(point: np.array):
        return -8.0

    expression = [
        DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
        SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
        SingleLayerVolumeTerm(Laplace2DKernel(), domain, heat_source_function, -1),
    ]

    problem = Problem(ProblemSolverType.BEM, expression, domain, analytical_solution)
    return problem


def init_poisson_neumann_bem():
    print("BEM for Neumann and Dirichlet mixed problem for Poisson equation...")

    print("Define boundary conditions...")

    def analytical_solution(point: np.array):
        return 2 * (point[0] - 0.5) ** 2 + 2 * (point[1] - 0.5) ** 2

    def dirichlet_boundary_value(point: np.array):
        return analytical_solution(point)

    def neumann_boundary_right_value(point: np.array):
        return 4 * point[0] - 2

    def neumann_boundary_left_value(point: np.array):
        return 2 - 4 * point[0]

    domain = SquareDomain2D(
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0]),
        [BoundaryConditionType.DIRICHLET, BoundaryConditionType.NEUMANN] * 2,
        [
            dirichlet_boundary_value,
            neumann_boundary_right_value,
            dirichlet_boundary_value,
            neumann_boundary_left_value,
        ],
        GlobalSettings.BORDER_ELEMENTS_COUNT,
    )

    print("Define heat source function...")

    def heat_source_function(point: np.array):
        return -8.0

    expression = [
        DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
        SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
        SingleLayerVolumeTerm(Laplace2DKernel(), domain, heat_source_function, -1),
    ]

    problem = Problem(ProblemSolverType.BEM, expression, domain, analytical_solution)
    return problem


def init_poisson_robin_bem():
    print("BEM for Robin problem for Poisson equation...")

    print("Define boundary conditions...")

    def analytical_solution(point: np.array):
        return point[0] + point[1] + 1.0

    def dirichlet_boundary_value(point: np.array):
        return analytical_solution(point)

    def robin_boundary_right_value(point: np.array):
        return (1.0 / analytical_solution(point), 0.0)

    def robin_boundary_left_value(point: np.array):
        return (-1.0 / analytical_solution(point), 0.0)

    domain = SquareDomain2D(
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0]),
        [BoundaryConditionType.DIRICHLET, BoundaryConditionType.ROBIN] * 2,
        [
            dirichlet_boundary_value,
            robin_boundary_right_value,
            dirichlet_boundary_value,
            robin_boundary_left_value,
        ],
        GlobalSettings.BORDER_ELEMENTS_COUNT,
    )

    print("Define heat source function...")

    def heat_source_function(point: np.array):
        return 0.0

    expression = [
        DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
        SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
        SingleLayerVolumeTerm(Laplace2DKernel(), domain, heat_source_function, -1),
    ]

    problem = Problem(ProblemSolverType.BEM, expression, domain, analytical_solution)
    return problem


def init_poisson_dirichlet_single_inclusion_bem():
    print("BEM for Dirichlet problem for Poisson equation with single inclusion...")

    print("Define boundary conditions...")

    K_TISSUE = 0.19  # W/m/^C
    K_MAX_TUMOR = 0.495  # W/m/^C
    INCLUSION_SIZE = 0.2
    INCLUSION_CENTER_X = 0.5
    INCLUSION_CENTER_Y = 0.5
    INCLUSION_RADIUS = INCLUSION_SIZE / 2.0 * np.sqrt(2.0)

    def boundary_value(point: np.array):
        return 2.0 * point[1]

    inclusion_domain = SquareDomain2D(
        np.array([INCLUSION_CENTER_X - INCLUSION_SIZE / 2.0, INCLUSION_CENTER_Y - INCLUSION_SIZE / 2.0]),
        np.array([INCLUSION_CENTER_X + INCLUSION_SIZE / 2.0, INCLUSION_CENTER_Y + INCLUSION_SIZE / 2.0]),
        [BoundaryConditionType.INCLUSION] * 4,
        [Utils.constant_one()] * 4,
        GlobalSettings.BORDER_ELEMENTS_COUNT,
    )

    domain = SquareDomain2D(
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0]),
        [BoundaryConditionType.DIRICHLET, BoundaryConditionType.NEUMANN] * 2,
        [boundary_value, Utils.constant_value(0.0)] * 2,
        GlobalSettings.BORDER_ELEMENTS_COUNT,
        [inclusion_domain],
    )

    def thermal_conductivity(point: np.array):
        result = K_TISSUE
        if inclusion_domain.is_point_inside_domain(point):
            result = (K_MAX_TUMOR - K_TISSUE) * np.cos(
                (0.5 * np.pi / INCLUSION_RADIUS**2)
                * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
            ) + K_TISSUE
        return result

    def thermal_conductivity_gradient(point: np.array):
        result = np.array([0.0, 0.0])
        if inclusion_domain.is_point_inside_domain(point):
            result[0] = (
                -(K_MAX_TUMOR - K_TISSUE)
                * np.sin(
                    (0.5 * np.pi / INCLUSION_RADIUS**2)
                    * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
                )
                * (np.pi / INCLUSION_RADIUS**2)
                * (point[0] - INCLUSION_CENTER_X)
            )
            result[1] = (
                -(K_MAX_TUMOR - K_TISSUE)
                * np.sin(
                    (0.5 * np.pi / INCLUSION_RADIUS**2)
                    * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
                )
                * (np.pi / INCLUSION_RADIUS**2)
                * (point[1] - INCLUSION_CENTER_Y)
            )
        return result

    def thermal_conductivity_laplacian(point: np.array):
        result = 0.0
        if inclusion_domain.is_point_inside_domain(point):
            result = ((K_TISSUE - K_MAX_TUMOR) * np.pi / INCLUSION_RADIUS**2) * (
                np.cos(
                    (0.5 * np.pi / INCLUSION_RADIUS**2)
                    * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
                )
                * (np.pi / INCLUSION_RADIUS**2 * (point[0] - INCLUSION_CENTER_X) ** 2)
                + np.sin(
                    (0.5 * np.pi / INCLUSION_RADIUS**2)
                    * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
                )
                + np.cos(
                    (0.5 * np.pi / INCLUSION_RADIUS**2)
                    * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
                )
                * (np.pi / INCLUSION_RADIUS**2 * (point[1] - INCLUSION_CENTER_X) ** 2)
                + np.sin(
                    (0.5 * np.pi / INCLUSION_RADIUS**2)
                    * ((point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2)
                )
            )
        return result

    print("Define heat source function...")

    expression = [
        DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
        SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
        SingleLayerInclusionTerm(
            Laplace2DKernel(), inclusion_domain, thermal_conductivity_gradient, thermal_conductivity_laplacian
        ),
    ]

    problem = Problem(ProblemSolverType.BEM, expression, domain)
    return problem


def init_poisson_dirichlet_cobem():
    print("CoBEM for Dirichlet problem for Poisson equation...")

    print("Define boundary conditions...")

    def analytical_solution(point: np.array):
        return 2 * (point[0] - 0.5) ** 2 + 2 * (point[1] - 0.5) ** 2

    def boundary_value(point: np.array):
        return analytical_solution(point)

    domain = SquareDomain2D(
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0]),
        [BoundaryConditionType.DIRICHLET] * 4,
        [boundary_value] * 4,
        GlobalSettings.BORDER_ELEMENTS_COUNT,
    )

    print("Define heat source function...")

    def heat_source_function(point: np.array):
        return -8.0

    expression = [
        SingleLayerCoBoundaryTerm(Laplace2DKernel(), domain),
        SingleLayerVolumeTerm(Laplace2DKernel(), domain, heat_source_function, -1),
    ]

    problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
    return problem


def init_poisson_neumann_cobem():
    print("CoBEM for Neumann and Dirichlet mixed problem for Poisson equation...")

    print("Define boundary conditions...")

    def analytical_solution(point: np.array):
        return 2.0 * (point[0] - 0.5) ** 2 + 2.0 * (point[1] - 0.5) ** 2

    def dirichlet_boundary_value(point: np.array):
        return analytical_solution(point)

    def neumann_boundary_right_value(point: np.array):
        return 4.0 * point[0] - 2.0

    def neumann_boundary_left_value(point: np.array):
        return 2.0 - 4.0 * point[0]

    domain = SquareDomain2D(
        np.array([0.0, 0.0]),
        np.array([1.0, 1.0]),
        [BoundaryConditionType.DIRICHLET, BoundaryConditionType.NEUMANN] * 2,
        [
            dirichlet_boundary_value,
            neumann_boundary_right_value,
            dirichlet_boundary_value,
            neumann_boundary_left_value,
        ],
        GlobalSettings.BORDER_ELEMENTS_COUNT,
    )

    print("Define heat source function...")

    def heat_source_function(point: np.array):
        return -8.0

    expression = [
        SingleLayerCoBoundaryTerm(Laplace2DKernel(), domain),
        SingleLayerVolumeCoBEMTerm(Laplace2DKernel(), domain, heat_source_function, -1),
    ]

    problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
    return problem


if "__main__" == __name__:
    problem = None
    if 1 == GlobalSettings.EXAMPLE_TYPE:
        problem = init_poisson_dirichlet_bem()
    elif 2 == GlobalSettings.EXAMPLE_TYPE:
        problem = init_poisson_neumann_bem()
    elif 3 == GlobalSettings.EXAMPLE_TYPE:
        problem = init_poisson_robin_bem()
    elif 4 == GlobalSettings.EXAMPLE_TYPE:
        problem = init_poisson_dirichlet_single_inclusion_bem()
    elif 5 == GlobalSettings.EXAMPLE_TYPE:
        problem = init_poisson_dirichlet_cobem()
    elif 6 == GlobalSettings.EXAMPLE_TYPE:
        problem = init_poisson_neumann_cobem()

    assert problem is not None

    problem.calculate()
    problem.plot()
