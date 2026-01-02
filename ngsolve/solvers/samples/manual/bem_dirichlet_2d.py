#!/usr/bin/env python3

from __future__ import annotations
import os
from enum import Enum
import math
from collections.abc import Callable
import numpy as np
import scipy
import csv
import statistics
import subprocess
import tempfile
import matplotlib.pyplot as plt
from matplotlib import cm, path, patches, transforms
from matplotlib.collections import PolyCollection
import meshio
from pathlib import Path


class GlobalSettings(object):
    CHART_STEPS = 25
    CHART_SKIP_BORDER_WIDTH = 0.0
    # TODO: Remove
    BORDER_ELEMENTS_COUNT = 10
    # TODO: Remove
    INCLUSION_ELEMENTS_COUNT = 5
    # TODO: Remove
    COBORDER_DEPTH = 1.2
    BORDER_ELEMENT_MAX_SIZE = 0.5
    INCLUSION_ELEMENTS_MAX_SIZE = 0.2
    COBORDER_SCALE = 2.0
    PLOT_ERROR = False
    PLOT_ISOLINES_COUNT = 10
    PLOT_ERROR_CONTOURS = False
    PLOT_DETAILS = False
    ERROR_TO_CSV = False
    ERROR_CSV_STEPS = 5
    """
        1 - Dirichlet BEM, 2 - Neumann BEM, 3 - Robin BEM, 4 - Single Inclusion Dirichlet BEM
        5 - Dirichlet CoBEM, 6 - Neumann CoBEM, 7 - Robin CoBEM, 8 - Single Inclusion Dirichlet CoBEM
        9 - Pennes Dirichlet BEM, 10 - Pennes Neumann BEM, 11 - Pennes Dirichlet CoBEM, 12 - Pennes Neumann CoBEM  
        13 - Pennes Dirichlet Hexagon CoBEM, 14 - Laplace Dirichlet Hexagon CoBEM, 
        15 - Laplace Dirichlet Convex BEM, 16 - Laplace Dirichlet Convex CoBEM,
        17 - Poisson Dirichlet Convex BEM, 18 - Poisson Dirichlet Convex CoBEM,
        19 - Laplace Dirichlet Single Inclusion Circular Convex BEM
    """
    EXAMPLE_TYPE = 19


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


class Constants:
    TWO_DIM = 2
    PLANE_SQUARE_DIM = 4


class Point2DInfo:
    def __init__(self):
        self.point = np.zeros(Constants.TWO_DIM)
        self.type = BoundaryConditionType.UNKNOWN
        self.normal = np.zeros(Constants.TWO_DIM)
        self.element = np.empty(0)
        self.value = 0.0
        self.robin_coeff = 0.0  # Robin condition is represented as 1 * q = robin_coeff * u + value
        self.is_right_corner = False


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
    def is_the_same_vector(vector1: np.array, vector2: np.array, eps: float):
        norm_vect1 = Utils.normalize(vector1)
        norm_vect2 = Utils.normalize(vector2)
        norm = np.linalg.norm(norm_vect1 - norm_vect2)
        result = norm < eps
        return result

    @staticmethod
    def normalize(vector: np.array):
        norm = np.linalg.norm(vector)
        if norm > 0.0 or norm < 0.0:
            vector /= norm
        return vector

    @staticmethod
    def get_normalized_vector(start_point: np.array, end_point: np.array):
        vector = end_point - start_point
        result = Utils.normalize(vector)
        return result

    @staticmethod
    def cross_product(vector1: np.array, vector2: np.array):
        if Constants.TWO_DIM == len(vector1) and Constants.TWO_DIM == len(vector2):
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

    # Calculate initial sign (expect all same sign for convex) based on cross product
    # (x2-x1)*(py-y1) - (y2-y1)*(px-x1) gives Z-component of cross product
    # Sign tells if point is left/right/on edge
    @staticmethod
    def is_point_within_convex_quadrilateral(point: np.array, vertices: np.array):
        result = True
        px, py = point
        num_vertices = len(vertices)
        p1 = vertices[0]
        p2 = vertices[1]
        initial_sign = np.sign((p2[0] - p1[0]) * (py - p1[1]) - (p2[1] - p1[1]) * (px - p1[0]))

        for i in range(1, num_vertices):
            p1 = vertices[i]
            p2 = vertices[(i + 1) % num_vertices]
            current_sign = np.sign((p2[0] - p1[0]) * (py - p1[1]) - (p2[1] - p1[1]) * (px - p1[0]))
            # If signs differ (and not zero), point is outside
            if current_sign != initial_sign and current_sign != 0:
                result = False
                break
        return result

    @staticmethod
    def is_square_side(begin: np.array, end: np.array, mesh_item: np.array, eps: float):
        assert Constants.TWO_DIM == len(begin)
        assert Constants.TWO_DIM == len(end)
        assert Constants.PLANE_SQUARE_DIM == len(mesh_item)

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

    @staticmethod
    def is_rhombus(points: np.array, eps: float):
        assert Constants.PLANE_SQUARE_DIM == len(points)

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

    @staticmethod
    def plot(x_limits: np.array, y_limits: np.array, functions: list, domain: Domain2D = None):
        assert Constants.TWO_DIM == len(x_limits)
        assert Constants.TWO_DIM == len(y_limits)
        assert len(functions) > 0

        functions_count = len(functions)
        fig = plt.figure(figsize=plt.figaspect(0.5))
        x_data_linear = np.linspace(x_limits[0], x_limits[1], GlobalSettings.CHART_STEPS)
        y_data_linear = np.linspace(y_limits[0], y_limits[1], GlobalSettings.CHART_STEPS)
        x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)

        for index, function in enumerate(functions):
            axis = fig.add_subplot(1, functions_count, index + 1, projection="3d")
            z_data = np.empty_like(x_data)
            z_data.fill(np.nan)
            for i in range(x_data.shape[0]):
                for j in range(x_data.shape[1]):
                    simple_point = (x_data[i, j], y_data[i, j])
                    point = np.array(simple_point)
                    if domain is not None and domain.contains_point(simple_point) or domain is None:
                        z_data[i, j] = function(point)
            surf = axis.plot_surface(x_data, y_data, z_data, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            fig.colorbar(surf, shrink=0.5, aspect=10)

        plt.show()

    @staticmethod
    def plot_contour(x_limits: np.array, y_limits: np.array, function, domain: Domain2D = None):
        assert Constants.TWO_DIM == len(x_limits)
        assert Constants.TWO_DIM == len(y_limits)
        assert function is not None

        fig, ax = plt.subplots()
        x_data_linear = np.linspace(x_limits[0], x_limits[1], GlobalSettings.CHART_STEPS)
        y_data_linear = np.linspace(y_limits[0], y_limits[1], GlobalSettings.CHART_STEPS)
        x_data, y_data = np.meshgrid(x_data_linear, y_data_linear)
        z_data = np.empty_like(x_data)
        z_data.fill(np.nan)
        for i in range(x_data.shape[0]):
            for j in range(x_data.shape[1]):
                simple_point = (x_data[i, j], y_data[i, j])
                point = np.array(simple_point)
                if domain is not None and domain.contains_point(simple_point) or domain is None:
                    z_data[i, j] = function(point)
        cs = ax.contour(x_data, y_data, z_data, levels=GlobalSettings.PLOT_ISOLINES_COUNT)
        ax.clabel(cs, inline=True, fontsize=8)

        plt.show()

    @staticmethod
    def plot_2d_domain(domain: Domain2D):
        assert domain is not None
        fig, ax = plt.subplots()
        patch = patches.PathPatch(domain.get_matplot_coborder(), facecolor="lightgreen", edgecolor="green", lw=2)
        ax.add_patch(patch)
        patch = patches.PathPatch(domain.get_matplot_border(), facecolor="lightblue", edgecolor="blue", lw=2)
        ax.add_patch(patch)
        for subdomain in domain.get_subdomains():
            patch = patches.PathPatch(
                subdomain.get_matplot_border(), facecolor="lightyellow", edgecolor="yellow", lw=2
            )
            ax.add_patch(patch)
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        plt.show()

    @staticmethod
    def plot_2d_mesh(mesh_info: list):
        assert mesh_info is not None
        mesh = [point_info.element for point_info in mesh_info]
        min_x = min(vertex[0] for quad in mesh for vertex in quad)
        min_y = min(vertex[1] for quad in mesh for vertex in quad)
        max_x = max(vertex[0] for quad in mesh for vertex in quad)
        max_y = max(vertex[1] for quad in mesh for vertex in quad)
        x_border = (max_x - min_x) * 0.5
        y_border = (max_y - min_y) * 0.5
        fig, ax = plt.subplots()
        col = PolyCollection(mesh, closed=True, facecolors="none")
        ax.add_collection(col)
        ax.set_xlim(min_x - x_border, max_x + x_border)
        ax.set_ylim(min_y - y_border, max_y + y_border)
        plt.show()

    @staticmethod
    def get_right_side_normal(begin: np.array, end: np.array):
        vector = Utils.get_normalized_vector(begin, end)
        result = np.array([vector[1], -1.0 * vector[0]])
        return result

    @staticmethod
    def get_next_index(current: int, list: np.array):
        assert list is not None
        result = current + 1
        if len(list) == result:
            result = 0
        return result

    @staticmethod
    def get_line_chunk_point(index: int, count: int, beginning: np.array, end: np.array):
        assert count > 0
        assert index >= 0 and index <= count
        assert beginning is not None
        assert end is not None

        result = beginning + (end - beginning) * float(index) / float(count)
        return result

    @staticmethod
    def order_quad_ccw(points: list):
        pts = np.asarray(points, dtype=float)
        assert pts.shape == (Constants.PLANE_SQUARE_DIM, Constants.TWO_DIM), f"Actual shape '{pts.shape}'."

        # Compute centroid
        cx, cy = pts.mean(axis=0)

        # Compute angle of each point around centroid
        angles = np.arctan2(pts[:, 1] - cy, pts[:, 0] - cx)

        # Sort points by angle (anticlockwise)
        order = np.argsort(angles)

        return pts[order]

    @staticmethod
    def get_file_in_current_directory(relative_path: str):
        assert path is not None
        module_dir = Path(__file__).resolve().parent
        return (module_dir / relative_path).resolve()


class MeshLoader:
    TOP = "top"
    LEFT = "left"
    RIGHT = "right"
    BOTTOM = "bottom"

    BODY = "body"

    ORDER = [BOTTOM, RIGHT, TOP, LEFT]

    UNIT_SQUARE_FILE = "unit_square.geo"
    UNIT_CIRCLE_FILE = "unit_circle.geo"

    @staticmethod
    def generate_mesh_from_file(filename: str, max_element_size: float, scale: float = 1.0):
        geo_folder = Utils.get_file_in_current_directory("geometries")
        geo_file = (geo_folder / filename).resolve()
        assert os.path.isfile(geo_file), f"File '{geo_file}' does not exist"
        result = None
        with tempfile.NamedTemporaryFile(suffix=".msh") as tmp:
            process_result = subprocess.run(
                [
                    "gmsh",
                    geo_file,
                    "-2",
                    "-clmax",
                    str(max_element_size),
                    "-scale",
                    str(scale),
                    "-setnumber",
                    "Mesh.Algorithm",
                    "6",  # TODO: Add check if only quad mesh was generated
                    "-setnumber",
                    "Mesh.RecombineAll",
                    "1",
                    "-setnumber",
                    "Mesh.RecombinationAlgorithm",
                    "1",
                    "-o",
                    tmp.name,
                ],
                capture_output=True,
                text=True,
            )
            assert 0 == process_result.returncode, f"Mesh generation error: {process_result.stderr}"
            result = meshio.read(tmp.name)
        return result

    @staticmethod
    def order_conditions(conditions: dict):
        assert conditions is not None
        result = {k: conditions[k] for k in MeshLoader.ORDER if k in conditions}
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

    def is_corner_point(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def is_point_inside_domain(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

    def contains_point(self, point: np.array):
        result = False
        if self.is_point_inside_domain(point) or self.is_point_on_border(point):
            result = True
        return result

    def get_subdomains(self):
        return self.__subdomains

    def print_stats(self):
        raise NotImplementedError("Call to abstract method")


class Domain2D(Domain):
    POINT_LOCATION_EPSILON = 0.001

    def __init__(self, subdomains: list = []):
        super().__init__(subdomains)
        self.__border_polygon = None
        self.__coborder_polygon = None

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

    def _init_polygon(self, polygon: np.array):
        assert polygon is not None
        border_list = [(point[0], point[1]) for point in polygon]
        item = path.Path(border_list)
        self._set_polygon(item)

    def _set_polygon(self, polygon: path.Path):
        assert polygon is not None
        self.__border_polygon = polygon

    def _init_coborder_polygon(self, center: np.array):
        assert center is not None
        assert Constants.TWO_DIM == len(center)
        assert self.__border_polygon is not None
        center_transform = transforms.Affine2D().translate(-center[0], -center[1])
        coboundary_scale = transforms.Affine2D().scale(GlobalSettings.COBORDER_SCALE)
        coborder_item = (
            self.__border_polygon.transformed(center_transform)
            .transformed(coboundary_scale)
            .transformed(center_transform.inverted())
        )
        self._set_coborder_polygon(coborder_item)

    def _set_coborder_polygon(self, polygon: path.Path):
        assert polygon is not None
        self.__coborder_polygon = polygon

    def get_matplot_border(self):
        assert self.__border_polygon is not None
        return self.__border_polygon

    def get_matplot_coborder(self):
        assert self.__coborder_polygon is not None
        return self.__coborder_polygon

    @staticmethod
    def _get_unique_vertices(vertices):
        assert vertices is not None
        outcome = []
        for index, item in enumerate(vertices):
            if index < len(vertices) - 1:
                next_item = np.array(vertices[index + 1])
            else:
                next_item = np.array(vertices[0])
            current = np.array(item)
            if not Utils.is_the_same_point(current, next_item, Domain2D.POINT_LOCATION_EPSILON):
                outcome.append(current)
        result = np.array(outcome)
        return result

    def is_point_on_border(self, point: np.array):
        for border_point in self.get_border():
            if Utils.is_point_within_segment(
                point, border_point.element[0], border_point.element[1], Domain2D.POINT_LOCATION_EPSILON
            ):
                return True
        return False

    def is_point_inside_domain(self, point: np.array):
        assert self.__border_polygon is not None
        return self.__border_polygon.contains_point(tuple(point))

    def is_corner_point(self, point: np.array):
        result = False
        for item in self.get_border():
            if item.is_right_corner and Utils.is_the_same_point(
                item.element[1], point, PlainDomain2D.POINT_LOCATION_EPSILON
            ):
                result = True
                break
        return result


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
        self._init_polygon(self.__square)
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

    def is_point_inside_domain(self, point: np.array):
        result = Utils.is_point_within_square(point, self.__square, Domain2D.POINT_LOCATION_EPSILON)
        return result

    def is_corner_point(self, point: np.array):
        result = False
        for corner in self.__square:
            if Utils.distance(point, corner) < Domain2D.POINT_LOCATION_EPSILON:
                result = True
                break
        return result


# TODO: Remove not used domain
class HexagonalDomain2D(Domain2D):
    def __init__(
        self,
        bottom_left_point: np.array,
        top_right_point: np.array,
        width: float,
        conditions: list[BoundaryConditionType],
        values: list[Callable],
        side_elements_count: int,
        subdomains: list = [],
    ):
        super().__init__(subdomains)

        assert bottom_left_point is not None
        assert top_right_point is not None
        assert width > 0.0
        assert 6 == len(conditions)
        assert 6 == len(values)
        assert side_elements_count > 0

        side_points_count = side_elements_count + 1
        center_point = (bottom_left_point + top_right_point) / 2.0
        self.__border = np.empty(0)

        self.__edge_points = []
        self.__edge_points.append(bottom_left_point)
        self.__edge_points.append(np.array([top_right_point[0], bottom_left_point[1]]))
        self.__edge_points.append(np.array([center_point[0] + width / 2.0, center_point[1]]))
        self.__edge_points.append(top_right_point)
        self.__edge_points.append(np.array([bottom_left_point[0], top_right_point[1]]))
        self.__edge_points.append(np.array([center_point[0] - width / 2.0, center_point[1]]))

        for index in range(0, len(self.__edge_points)):
            begin = self.__edge_points[index]
            if index < len(self.__edge_points) - 1:
                end = self.__edge_points[index + 1]
            else:
                end = self.__edge_points[0]
            line_points = np.linspace(begin, end, side_points_count)
            self.__border = np.append(
                self.__border, np.array(Domain2D._process_points(line_points, conditions[index], values[index]))
            )

        self._init_polygon(self.__edge_points)

        # TODO: Add Mesh
        self.__mesh = []
        self.__mesh = np.array(self.__mesh)
        # TODO: Add square
        # self.__square = np.array([bottom_left_point, bottom_right_point, top_right_point, top_left_point])
        self._fill_coborder_elements()

    def _fill_coborder_elements(self):
        assert self.__border is not None
        assert self.__edge_points is not None

        def get_growth_vector(prev_start, corner, next_end, normal1, normal2):
            direction_previous_forward = Utils.get_normalized_vector(prev_start, corner)
            direction_next_backward = Utils.get_normalized_vector(next_end, corner)
            cross_prod_norms = Utils.cross_product(normal1, normal2)
            result = direction_previous_forward + direction_next_backward
            result = result * GlobalSettings.COBORDER_DEPTH / np.linalg.norm(cross_prod_norms)
            return result

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
            for corner in self.__edge_points:
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
                    corner_grow_vector = get_growth_vector(
                        self.__border[prev_index].element[0],
                        point_info.element[0],
                        point_info.element[1],
                        self.__border[prev_index].normal,
                        point_info.normal,
                    )
                    elements.append(point_info.element[0] + corner_grow_vector)
                    if 1 == len(corner_points):
                        elements.append(point_info.element[1] + point_info.normal * GlobalSettings.COBORDER_DEPTH)

                if Utils.is_the_same_point(elements[1], corner_point, Domain2D.POINT_LOCATION_EPSILON) or 2 == len(
                    corner_points
                ):
                    next_index = index + 1
                    if next_index >= len(self.__border):
                        next_index = 0
                    corner_grow_vector = get_growth_vector(
                        point_info.element[0],
                        point_info.element[1],
                        self.__border[next_index].element[1],
                        point_info.normal,
                        self.__border[next_index].normal,
                    )
                    elements.append(point_info.element[1] + corner_grow_vector)
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

    def get_edge_points(self):
        return self.__edge_points


# GUI usage tutorial
# https://www.youtube.com/watch?v=kk5DHIZa21k&list=PLLaFJ14_gbrM9nlZfYcuXF40QxZb-ygfY
class GmshDomain2D(Domain2D):
    POINT_LOCATION_EPSILON = 0.001
    GMSH_CELL_TYPE_LINE = "line"
    GMSH_CELL_TYPE_QUAD = "quad"

    def __init__(self, loaded_mesh, conditions: dict, center: np.array = None, subdomains: list = []):
        super().__init__(subdomains)
        assert loaded_mesh is not None
        assert conditions is not None
        assert center is None or Constants.TWO_DIM == len(center)
        self.__border = []
        self.__coborder = []
        self.__mesh = []
        if center is None:
            self.__center = np.zeros(Constants.TWO_DIM)
        else:
            self.__center = center
        self.__init_border(conditions, loaded_mesh)
        self._init_coborder_polygon(self.__center)
        self.__process_coborder_points()
        self.__process_mesh(loaded_mesh)

    def __init_border(self, conditions, mesh):
        polygon = []
        for physical_name, (condition, value_function) in conditions.items():
            assert physical_name in mesh.field_data, f"Physical group '{physical_name}' not found in mesh."
            physical_tag, physical_dim = mesh.field_data[physical_name]
            for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data.get("gmsh:physical", [])):
                if cell_block.dim != physical_dim or cell_block.type != GmshDomain2D.GMSH_CELL_TYPE_LINE:
                    continue
                mask = phys_tags == physical_tag
                selected = cell_block.data[mask]
                for item in selected:
                    start_point = np.resize(mesh.points[item[0]], (Constants.TWO_DIM,)) + self.__center
                    end_point = np.resize(mesh.points[item[1]], (Constants.TWO_DIM,)) + self.__center
                    if 0 == len(polygon):
                        polygon.append(start_point)
                    polygon.append(end_point)
                    self.__border = np.append(
                        self.__border, Domain2D._process_points([start_point, end_point], condition, value_function)
                    )
        self.__init_corner_points()
        self._init_polygon(polygon)

    # TODO: Remove duplicate from PlainDomain2D
    def __init_corner_points(self):
        assert len(self.__border) > 0, "Border points should be initialized"
        for index, element in enumerate(self.__border):
            if index < len(self.__border) - 1:
                next_element = self.__border[index + 1]
            else:
                next_element = self.__border[0]
            if not Utils.is_the_same_vector(element.normal, next_element.normal, PlainDomain2D.POINT_LOCATION_EPSILON):
                element.is_right_corner = True

    # TODO: Remove duplicate from PlainDomain2D
    def __process_coborder_points(self):
        assert self.__border is not None
        border_vertices = Domain2D._get_unique_vertices(self.get_matplot_border().vertices)
        coborder_vertices = Domain2D._get_unique_vertices(self.get_matplot_coborder().vertices)
        for index in range(len(self.__border)):
            point_info = self.__border[index]
            item = Point2DInfo()
            item.point = point_info.point
            item.type = point_info.type
            item.normal = point_info.normal
            item.value = point_info.value
            item.robin_coeff = point_info.robin_coeff

            next_index = index + 1
            if next_index == len(self.__border):
                next_index = 0

            elements = []
            elements.append(coborder_vertices[index])
            elements.append(coborder_vertices[next_index])
            elements.append(border_vertices[next_index])
            elements.append(border_vertices[index])
            item.element = elements
            self.__coborder.append(item)

        self.__coborder = np.array(self.__coborder)

    def __process_mesh(self, mesh):
        assert mesh is not None
        for cell_block in mesh.cells:
            if GmshDomain2D.GMSH_CELL_TYPE_QUAD == cell_block.type:
                for cell in cell_block.data:
                    quad = []
                    for index in cell:
                        quad_vertex = np.resize(mesh.points[index], (Constants.TWO_DIM,)) + self.__center
                        quad.append(quad_vertex)
                    quad = Utils.order_quad_ccw(quad)
                    point_info = Point2DInfo()
                    point_info.point = 0.5 * (quad[0] + quad[2])
                    point_info.type = BoundaryConditionType.INCLUSION
                    point_info.value = 0.0
                    point_info.element = quad
                    self.__mesh.append(point_info)
        self.__mesh = np.array(self.__mesh)

    def get_border(self):
        return self.__border

    def get_mesh(self):
        return self.__mesh

    def get_coborder(self):
        return self.__coborder

    # TODO: Remove duplicate from PlainDomain2D
    def print_stats(self):
        # TODO: Add max side size and min side size for all meshes and coborder width
        print("---")
        print("GmshDomain2D stats:")
        print("---")
        border_elements_lengths = [Utils.distance(item.element[0], item.element[1]) for item in self.__border]
        print(
            f"""Border elements 
              count: {len(self.__border)}
              min length: {min(border_elements_lengths)}
              median length: {statistics.median(border_elements_lengths)}
              max length: {max(border_elements_lengths)}"""
        )
        mesh_elements_lengths = [
            Utils.distance(start, end)
            for item in self.__mesh
            for start, end in zip(item.element, np.vstack([item.element[1:], item.element[:1]]))
        ]
        print(
            f"""Quad elements 
              count: {len(self.__mesh)}
              min side length: {min(mesh_elements_lengths)}
              median side length: {statistics.median(mesh_elements_lengths)}
              max side length: {max(mesh_elements_lengths)}"""
        )
        print("---")
        for domain in self.get_subdomains():
            domain.print_stats()


class PlainDomain2D(Domain2D):
    POINT_LOCATION_EPSILON = 0.001

    def __init__(
        self, paths: list[path.Path], conditions: list, functions: list, center: np.array = None, subdomains: list = []
    ):
        super().__init__(subdomains)
        assert paths is not None
        assert center is None or Constants.TWO_DIM == len(center)
        assert len(paths) == len(conditions)
        assert len(paths) == len(functions)
        if center is None:
            self.__center = np.zeros(Constants.TWO_DIM)
        else:
            self.__center = center
        self.__border = []
        list = []
        for index, item in enumerate(paths):
            item.to_polygons(closed_only=False)
            verts = item.vertices
            max_dist = max(Utils.distance(verts[i], verts[i - 1]) for i in range(1, len(verts)))
            steps = math.ceil(max_dist / GlobalSettings.BORDER_ELEMENT_MAX_SIZE)
            new_item = item.interpolated(steps)
            list.append(new_item)
            points = self.__process_points(new_item, conditions[index], functions[index])
            self.__border = np.append(self.__border, points)
        self.__init_corner_points()
        item = path.Path.make_compound_path(*list)
        self._set_polygon(item)
        self._init_coborder_polygon(self.__center)
        self.__process_coborder_points()
        self.__process_mesh()

    def __init_corner_points(self):
        assert len(self.__border) > 0, "Border points should be initialized"
        for index, element in enumerate(self.__border):
            if index < len(self.__border) - 1:
                next_element = self.__border[index + 1]
            else:
                next_element = self.__border[0]
            if not Utils.is_the_same_vector(element.normal, next_element.normal, PlainDomain2D.POINT_LOCATION_EPSILON):
                element.is_right_corner = True

    @staticmethod
    def __process_points(path_element: path.Path, type: BoundaryConditionType, value_function: Callable):
        assert path_element is not None
        assert value_function is not None
        result = []
        verts = Domain2D._get_unique_vertices(path_element.vertices)
        for index in range(len(verts) - 1):
            item = Point2DInfo()
            item.point = (verts[index + 1] + verts[index]) / 2.0
            item.type = type
            border_vector = verts[index + 1] - verts[index]
            assert np.linalg.norm(border_vector) != 0, "Points of border element should be different"
            border_vector /= np.linalg.norm(border_vector)

            # TODO: Check if it is always external normal
            item.normal = np.array([border_vector[1], -border_vector[0]])
            item.element = np.array([verts[index], verts[index + 1]])
            if BoundaryConditionType.ROBIN == type:
                item.robin_coeff, item.value = value_function(item.point)
            else:
                item.value = value_function(item.point)
            result.append(item)
        return result

    def __process_coborder_points(self):
        assert self.__border is not None
        self.__coborder = []
        border_vertices = Domain2D._get_unique_vertices(self.get_matplot_border().vertices)
        coborder_vertices = Domain2D._get_unique_vertices(self.get_matplot_coborder().vertices)
        for index in range(len(self.__border)):
            point_info = self.__border[index]
            item = Point2DInfo()
            item.point = point_info.point
            item.type = point_info.type
            item.normal = point_info.normal
            item.value = point_info.value
            item.robin_coeff = point_info.robin_coeff

            next_index = index + 1
            if next_index == len(self.__border):
                next_index = 0

            elements = []
            elements.append(coborder_vertices[index])
            elements.append(coborder_vertices[next_index])
            elements.append(border_vertices[next_index])
            elements.append(border_vertices[index])
            item.element = elements
            self.__coborder.append(item)

        self.__coborder = np.array(self.__coborder)

    # TODO: Implement for non-convex polygons
    def __process_mesh(self):
        assert self.__border is not None
        assert self.__center is not None
        self.__mesh = []
        for point in self.__border:
            mesh = PlainDomain2D.__quadragulate_triangle(point.element[0], self.__center, point.element[1])
            for quad in mesh:
                point_info = Point2DInfo()
                point_info.point = 0.5 * (quad[0] + quad[2])
                point_info.type = BoundaryConditionType.INCLUSION
                point_info.value = 0.0
                point_info.element = np.array(quad)
                self.__mesh.append(point_info)
        self.__mesh = np.array(self.__mesh)

    @staticmethod
    def __quadragulate_triangle(vertex1: np.array, vertex2: np.array, vertex3: np.array):
        result = []
        left_side_len = Utils.distance(vertex1, vertex2)
        right_side_len = Utils.distance(vertex3, vertex2)
        slices_count = math.ceil(max(left_side_len, right_side_len) / GlobalSettings.INCLUSION_ELEMENTS_MAX_SIZE)
        for slice_index in range(slices_count - 1):
            bottom_line_left = Utils.get_line_chunk_point(slice_index, slices_count, vertex1, vertex2)
            bottom_line_right = Utils.get_line_chunk_point(slice_index, slices_count, vertex3, vertex2)
            top_line_left = Utils.get_line_chunk_point(slice_index + 1, slices_count, vertex1, vertex2)
            top_line_right = Utils.get_line_chunk_point(slice_index + 1, slices_count, vertex3, vertex2)
            bottom_side_len = Utils.distance(bottom_line_left, bottom_line_right)
            chunks_count = math.ceil(bottom_side_len / GlobalSettings.INCLUSION_ELEMENTS_MAX_SIZE)
            for chunk_index in range(chunks_count):
                quadrilateral = [
                    Utils.get_line_chunk_point(chunk_index, chunks_count, bottom_line_left, bottom_line_right),
                    Utils.get_line_chunk_point(chunk_index + 1, chunks_count, bottom_line_left, bottom_line_right),
                    Utils.get_line_chunk_point(chunk_index + 1, chunks_count, top_line_left, top_line_right),
                    Utils.get_line_chunk_point(chunk_index, chunks_count, top_line_left, top_line_right),
                ]
                result.append(quadrilateral)

        # Process small triangle on top by splitting it using middle lines
        # finding point inside inner middle line based triangle and
        # joining middle point with vertices and medians
        vertex_left = Utils.get_line_chunk_point(slices_count - 1, slices_count, vertex1, vertex2)
        vertex_right = Utils.get_line_chunk_point(slices_count - 1, slices_count, vertex3, vertex2)
        mid_left_point = (vertex_left + vertex2) * 0.5
        mid_right_point = (vertex_right + vertex2) * 0.5
        mid_bottom_point = (vertex_left + vertex_right) * 0.5

        mid_mid_point_top = (mid_left_point + mid_right_point) * 0.5
        bulls_eye_point = (mid_bottom_point + mid_mid_point_top) * 0.5

        result.append([bulls_eye_point, mid_right_point, vertex2, mid_left_point])
        result.append([bulls_eye_point, mid_left_point, vertex_left, mid_bottom_point])
        result.append([bulls_eye_point, mid_bottom_point, vertex_right, mid_right_point])

        return result

    def get_border(self):
        return self.__border

    def get_mesh(self):
        return self.__mesh

    def get_coborder(self):
        return self.__coborder

    def print_stats(self):
        # TODO: Add max side size and min side size for all meshes
        print("PlainDomain2D stats:")
        print(f"Border elements count: {len(self.__border)}")
        print(f"Quads count: {len(self.__mesh)}")
        print("---")

        for domain in self.get_subdomains():
            domain.print_stats()


class KernelValueType(Enum):
    SCALAR = 1
    DX = 2
    DY = 3
    DXX = 4
    DXY = 5
    DYY = 6
    GRADIENT = 7
    GRADIENT_DX = 8
    GRADIENT_DY = 9


class Kernel:
    def value(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def dx(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def dy(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def dxx(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def dxy(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def dyy(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Call to abstract method")

    def grad(self, point_x: np.array, point_y: np.array):
        return np.array([self.dx(point_x, point_y), self.dy(point_x, point_y)])

    def grad_dx(self, point_x: np.array, point_y: np.array):
        return np.array([self.dxx(point_x, point_y), self.dxy(point_x, point_y)])

    def grad_dy(self, point_x: np.array, point_y: np.array):
        return np.array([self.dxy(point_x, point_y), self.dyy(point_x, point_y)])

    def value_of(self, type: KernelValueType, point_x: np.array, point_y: np.array):
        if KernelValueType.SCALAR == type:
            result = self.value(point_x, point_y)
        elif KernelValueType.DX == type:
            result = self.dx(point_x, point_y)
        elif KernelValueType.DY == type:
            result = self.dy(point_x, point_y)
        elif KernelValueType.DXX == type:
            result = self.dxx(point_x, point_y)
        elif KernelValueType.DXY == type:
            result = self.dxy(point_x, point_y)
        elif KernelValueType.DYY == type:
            result = self.dyy(point_x, point_y)
        elif KernelValueType.GRADIENT == type:
            result = self.grad(point_x, point_y)
        elif KernelValueType.GRADIENT_DX == type:
            result = self.grad_dx(point_x, point_y)
        elif KernelValueType.GRADIENT_DY == type:
            result = self.grad_dy(point_x, point_y)
        else:
            assert "Unknown kernel value type"
        return result


class Laplace2DKernel(Kernel):
    def value(self, point_x: np.array, point_y: np.array):
        return -0.25 / np.pi * np.log(Utils.squared_distance(point_x, point_y))

    def dx(self, point_x: np.array, point_y: np.array):
        result = -0.25 / np.pi / Utils.squared_distance(point_x, point_y) * -2.0 * (point_x[0] - point_y[0])
        return result

    def dy(self, point_x: np.array, point_y: np.array):
        result = -0.25 / np.pi / Utils.squared_distance(point_x, point_y) * -2.0 * (point_x[1] - point_y[1])
        return result

    def dxx(self, point_x: np.array, point_y: np.array):
        distance_forth = Utils.squared_distance(point_x, point_y) ** 2
        result = -0.5 / np.pi / distance_forth * ((point_x[0] - point_y[0]) ** 2 - (point_x[1] - point_y[1]) ** 2)
        return result

    def dxy(self, point_x: np.array, point_y: np.array):
        distance_forth = Utils.squared_distance(point_x, point_y) ** 2
        result = -1.0 / np.pi / distance_forth * (point_x[0] - point_y[0]) * (point_x[1] - point_y[1])
        return result

    def dyy(self, point_x: np.array, point_y: np.array):
        result = -1.0 * self.dxx(point_x, point_y)
        return result


# https://en.wikipedia.org/wiki/Green%27s_function#Table_of_Green's_functions
class Pennes2DKernel(Kernel):

    def __init__(self, k_square: float):
        super().__init__()
        assert k_square > 0
        self.__k = np.sqrt(k_square)

    def value(self, point_x: np.array, point_y: np.array):
        result = -0.5 / np.pi * scipy.special.kn(0, self.__k * Utils.distance(point_x, point_y))
        return result

    def dx(self, point_x: np.array, point_y: np.array):
        r = Utils.distance(point_x, point_y)
        result = (
            -0.5
            / np.pi
            * scipy.special.kvp(0, self.__k * r, 1)
            * 0.5
            * self.__k
            / r
            * -2.0
            * (point_x[0] - point_y[0])
        )
        return result

    def dy(self, point_x: np.array, point_y: np.array):
        r = Utils.distance(point_x, point_y)
        result = (
            -0.5
            / np.pi
            * scipy.special.kvp(0, self.__k * r, 1)
            * 0.5
            * self.__k
            / r
            * -2.0
            * (point_x[1] - point_y[1])
        )
        return result

    def dxx(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Implement me")

    def dxy(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Implement me")

    def dyy(self, point_x: np.array, point_y: np.array):
        raise NotImplementedError("Implement me")


class Integrator:
    def __init__(self, kernel: Kernel):
        assert kernel is not None
        self.__kernel = kernel

    @property
    def kernel(self):
        return self.__kernel

    # TODO: Deprecated. Use segment_of instead
    def segment(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        raise NotImplementedError("Call to abstract method")

    # TODO: Deprecated. Use segment_of instead
    def segment_grad(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        raise NotImplementedError("Call to abstract method")

    def segment_of(
        self, type: KernelValueType, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array
    ):
        raise NotImplementedError("Call to abstract method")

    # TODO: Deprecated. Use square_of instead
    def square(self, function: Callable, point: np.array, square: np.array):
        raise NotImplementedError("Call to abstract method")

    # TODO: Deprecated. Use square_of instead
    def square_grad(self, function: Callable, point: np.array, square: np.array):
        raise NotImplementedError("Call to abstract method")

    def square_of(self, type: KernelValueType, function: Callable, point: np.array, square: np.array):
        raise NotImplementedError("Call to abstract method")

    def trapezoid_of(self, type: KernelValueType, function: Callable, point: np.array, trapezoid: np.array):
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

    def segment_of(
        self, type: KernelValueType, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array
    ):
        result = 0.0
        real_nodes, real_weights = self._get_segment_nodes_and_weights(point, min_limit, max_limit)
        if type in [KernelValueType.GRADIENT, KernelValueType.GRADIENT_DX, KernelValueType.GRADIENT_DY]:
            for real_node, weight in zip(real_nodes, real_weights):
                result += weight * np.dot(self.kernel.grad(point, real_node), function(real_node))
        else:
            for real_node, weight in zip(real_nodes, real_weights):
                result += weight * self.kernel.value_of(type, point, real_node) * function(real_node)
        return result

    # TODO: Deprecated. Use segment_of instead
    def segment(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        result = self.segment_of(KernelValueType.SCALAR, function, point, min_limit, max_limit)
        return result

    # TODO: Deprecated. Use segment_of instead
    def segment_grad(self, function: Callable, point: np.array, min_limit: np.array, max_limit: np.array):
        result = self.segment_of(KernelValueType.GRADIENT, function, point, min_limit, max_limit)
        return result

    def _get_square_nodes_and_weights(self, point: np.array, square: np.array):
        count = self.__nominal_number_of_points
        if Utils.is_point_within_square(point, square, Integrator2D.EPSILON):
            count = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(count)
        (real_nodes, real_weights) = self._convert_leggauss_to_square(nodes, weights, square)
        return (real_nodes, real_weights)

    # TODO: Deprecated. Use square_of instead
    def square(self, function: Callable, point: np.array, square: np.array):
        result = self.square_of(KernelValueType.SCALAR, function, point, square)
        return result

    # TODO: Deprecated. Use square_of instead
    def square_grad(self, function: Callable, point: np.array, square: np.array):
        result = self.square_of(KernelValueType.GRADIENT, function, point, square)
        return result

    def square_of(self, type: KernelValueType, function: Callable, point: np.array, square: np.array):
        result = 0.0
        real_nodes, real_weights = self._get_square_nodes_and_weights(point, square)
        if type in [KernelValueType.GRADIENT, KernelValueType.GRADIENT_DX, KernelValueType.GRADIENT_DY]:
            for real_node, weight in zip(real_nodes, real_weights):
                result += weight * np.dot(self.kernel.grad(point, real_node), function(real_node))
        else:
            for real_node, weight in zip(real_nodes, real_weights):
                result += weight * self.kernel.value_of(type, point, real_node) * function(real_node)
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

    # TODO: Deprecated. Use trapezoid_of instead
    def trapezoid(self, function: Callable, point: np.array, trapezoid: np.array):
        result = self.trapezoid_of(KernelValueType.SCALAR, function, point, trapezoid)
        return result

    # TODO: Deprecated. Use trapezoid_of instead
    def trapezoid_grad(self, function: Callable, point: np.array, trapezoid: np.array):
        result = self.trapezoid_of(KernelValueType.GRADIENT, function, point, trapezoid)
        return result

    def trapezoid_of(self, type: KernelValueType, function: Callable, point: np.array, trapezoid: np.array):
        assert 4 == len(trapezoid)
        result = 0.0
        if Utils.is_rhombus(trapezoid, Integrator2D.EPSILON):
            result = self.square_of(type, function, point, trapezoid)
            return result
        nodes, weights, coeff = self._get_trapezoid_nodes_weights_and_coeff(point, trapezoid)
        for node, weight in zip(nodes, weights):
            real_node = np.zeros(Integrator2D.DIM)
            a, b, c, d, e, f = (coeff[1], coeff[2], coeff[3], coeff[5], coeff[6], coeff[7])
            const_x, const_y = (coeff[0], coeff[4])
            real_node[0] = const_x + a * node[0] + b * node[1] + c * node[0] * node[1]
            real_node[1] = const_y + d * node[0] + e * node[1] + f * node[0] * node[1]
            jacobian = (a * e - b * d) + (a * f - c * d) * node[0] + (c * e - b * f) * node[1]
            if type in [KernelValueType.GRADIENT, KernelValueType.GRADIENT_DX, KernelValueType.GRADIENT_DY]:
                result += weight * np.dot(self.kernel.value_of(type, point, real_node), function(real_node)) * jacobian
            else:
                result += weight * self.kernel.value_of(type, point, real_node) * function(real_node) * jacobian
        return result

    def convex_quadrilateral_of(self, type: KernelValueType, function: Callable, point: np.array, vertices: np.array):
        assert Constants.PLANE_SQUARE_DIM == len(vertices)
        result = 0.0
        unit_square = np.array(
            [np.array([-1.0, -1.0]), np.array([1.0, -1.0]), np.array([1.0, 1.0]), np.array([-1.0, 1.0])]
        )
        matrix = None
        right_side = np.empty(0)
        for unit, actual in zip(unit_square, vertices):
            matrix_row = np.array([1.0, unit[0], unit[1], unit[0] * unit[1], 0.0, 0.0, 0.0, 0.0])
            if matrix is None:
                matrix = matrix_row
            else:
                matrix = np.vstack((matrix, matrix_row))
            matrix_row = np.array([0.0, 0.0, 0.0, 0.0, 1.0, unit[0], unit[1], unit[0] * unit[1]])
            matrix = np.vstack((matrix, matrix_row))
            right_side = np.append(right_side, [actual[0], actual[1]])
        coeff = np.linalg.solve(matrix, right_side)
        a, b, c, d, e, f = (coeff[1], coeff[2], coeff[3], coeff[5], coeff[6], coeff[7])
        const_x, const_y = (coeff[0], coeff[4])

        interpolation_order = self.__nominal_number_of_points
        if Utils.is_point_within_convex_quadrilateral(point, vertices):
            interpolation_order = self.__singular_number_of_points
        nodes, weights = np.polynomial.legendre.leggauss(interpolation_order)

        for x1, wx in zip(nodes, weights):
            for y1, wy in zip(nodes, weights):
                node = np.array([const_x + a * x1 + b * y1 + c * x1 * y1, const_y + d * x1 + e * y1 + f * x1 * y1])
                jacobian = (a * e - b * d) + (a * f - c * d) * x1 + (c * e - b * f) * y1
                if type in [KernelValueType.GRADIENT, KernelValueType.GRADIENT_DX, KernelValueType.GRADIENT_DY]:
                    result += wx * wy * np.dot(self.kernel.value_of(type, point, node), function(node)) * jacobian
                else:
                    result += wx * wy * self.kernel.value_of(type, point, node) * function(node) * jacobian
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

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
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

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0
        for boundary_item in self.domain.get_border():
            if boundary_item.type in [BoundaryConditionType.DIRICHLET, BoundaryConditionType.ROBIN]:
                res = self.__integrator.segment_of(
                    KernelValueType.SCALAR,
                    Utils.constant_one(),
                    point_info.point,
                    boundary_item.element[0],
                    boundary_item.element[1],
                )
                coefficients = np.append(coefficients, [res])
            elif BoundaryConditionType.NEUMANN == boundary_item.type:
                right_side_ret += boundary_item.value * self.__integrator.segment_of(
                    KernelValueType.SCALAR,
                    Utils.constant_one(),
                    point_info.point,
                    boundary_item.element[0],
                    boundary_item.element[1],
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
            integral_value = self.__integrator.segment_of(
                KernelValueType.SCALAR, Utils.constant_one(), point, boundary_item.element[0], boundary_item.element[1]
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

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0
        for boundary_item in self.domain.get_border():
            is_same_point = Utils.is_the_same_point(point_info.point, boundary_item.point, DoubleLayerBoundaryTerm.EPS)
            if BoundaryConditionType.DIRICHLET == boundary_item.type:
                right_side_ret += boundary_item.value * self.__integrator.segment_of(
                    KernelValueType.GRADIENT,
                    Utils.constant_value(boundary_item.normal),
                    point_info.point,
                    boundary_item.element[0],
                    boundary_item.element[1],
                )
                if is_same_point:
                    right_side_ret += 0.5 * boundary_item.value
            elif boundary_item.type in [BoundaryConditionType.NEUMANN, BoundaryConditionType.ROBIN]:
                res = self.__integrator.segment_of(
                    KernelValueType.GRADIENT,
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
            integral_value = self.__integrator.segment_of(
                KernelValueType.GRADIENT,
                Utils.constant_value(boundary_item.normal),
                point,
                boundary_item.element[0],
                boundary_item.element[1],
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

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = self.value_of(KernelValueType.SCALAR, point_info.point)
        return (coefficients, right_side_ret)

    def calculate_for_robin(self, point_info: Point2DInfo):
        coefficients = np.empty(0)
        return (coefficients, 0.0)

    def propagate_solution(self, solution: np.array):
        return solution

    def value_of(self, type: KernelValueType, point: np.array):
        assert point is not None
        result = 0.0
        for point_info in self.domain.get_mesh():
            result += self.__integrator.convex_quadrilateral_of(type, self.__value_function, point, point_info.element)
        result *= self.sign
        return result

    def value(self, point: np.array):
        assert point is not None
        result = self.value_of(KernelValueType.SCALAR, point)
        return result


class SingleLayerVolumeCoBEMTerm(SingleLayerVolumeTerm):

    EPS = 0.001

    def __init__(self, kernel: Kernel, domain: Domain, value_function: Callable, sign: int = 1):
        super().__init__(kernel, domain, value_function, sign)

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
        assert point_info is not None
        coefficients = np.empty(0)
        # TODO: Add Robin condition
        if BoundaryConditionType.NEUMANN == point_info.type:
            right_side_ret = self.value_grad_of(KernelValueType.GRADIENT, point_info.point, point_info.normal)
        else:
            right_side_ret = super().value_of(KernelValueType.SCALAR, point_info.point)
        return (coefficients, right_side_ret)

    def value_grad_of(self, type: KernelValueType, point: np.array, norm: np.array):
        assert point is not None
        assert type in [KernelValueType.GRADIENT, KernelValueType.GRADIENT_DX, KernelValueType.DY]
        result = 0.0

        def normal_function(x: np.array):
            result = self._get_value_function()(x) * norm
            return result

        for point_info in self.domain.get_mesh():
            result += self._get_integrator().convex_quadrilateral_of(type, normal_function, point, point_info.element)
        result *= self.sign
        return result


class SingleLayerCoBEMTerm(ExpressionTerm):
    NOMINAL_INTEGRATION_POINTS = 4
    SINGULARITY_INTEGRATION_POINTS = 6
    EPS = 0.001

    def __init__(self, kernel: Kernel, domain: Domain, sign: int = 1):
        super().__init__(kernel, domain, sign)
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)
        self.__integrator = Integrator2D(
            kernel,
            SingleLayerCoBEMTerm.NOMINAL_INTEGRATION_POINTS,
            SingleLayerCoBEMTerm.SINGULARITY_INTEGRATION_POINTS,  # There should not be singularities
        )

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
        assert point_info is not None
        coefficients = np.empty(0)
        right_side_ret = 0.0
        for boundary_item in self.domain.get_coborder():
            is_same_point = Utils.is_the_same_point(point_info.point, boundary_item.point, SingleLayerCoBEMTerm.EPS)
            if point_info.type in [BoundaryConditionType.DIRICHLET, BoundaryConditionType.INCLUSION]:
                res = self.__integrator.trapezoid_of(
                    KernelValueType.SCALAR, Utils.constant_one(), point_info.point, boundary_item.element
                )
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value
            elif BoundaryConditionType.NEUMANN == point_info.type:
                res = self.__integrator.trapezoid_of(
                    KernelValueType.GRADIENT,
                    Utils.constant_value(boundary_item.normal),
                    point_info.point,
                    boundary_item.element,
                )
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value
            elif BoundaryConditionType.ROBIN == point_info.type:
                res = self.__integrator.trapezoid_of(
                    KernelValueType.GRADIENT,
                    Utils.constant_value(point_info.normal),
                    point_info.point,
                    boundary_item.element,
                )
                res -= point_info.robin_coeff * self.__integrator.trapezoid_of(
                    KernelValueType.SCALAR, Utils.constant_one(), point_info.point, boundary_item.element
                )
                coefficients = np.append(coefficients, [res])
                if is_same_point:
                    right_side_ret += boundary_item.value
            else:
                raise AttributeError("Not supported boundary element type")
        self.__unknown_count = len(coefficients)

        return (self.sign * coefficients, self.sign * right_side_ret)

    def calculate_for_robin(self, point_info: Point2DInfo):
        raise NotImplementedError("Should not be called")

    def propagate_solution(self, solution: np.array):
        self.__unknown_values = solution[: self.__unknown_count]
        return solution[self.__unknown_count :]

    def value(self, point: np.array):
        assert self.__unknown_count == len(self.__unknown_values), "Data should be calculated"
        result = 0.0
        unknown_index = 0

        for boundary_item in self.domain.get_coborder():
            integral_value = self.__integrator.trapezoid_of(
                KernelValueType.SCALAR, Utils.constant_one(), point, boundary_item.element
            )
            if boundary_item.type in [
                BoundaryConditionType.DIRICHLET,
                BoundaryConditionType.NEUMANN,
                BoundaryConditionType.ROBIN,
            ]:
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
    EPS = 0.001

    def __init__(
        self,
        kernel: Kernel,
        domain: Domain,
        coeff_function: Callable,
        grad_function: Callable,
        sign: int = 1,
    ):
        super().__init__(kernel, domain, sign)
        assert coeff_function is not None
        assert grad_function is not None

        def adjusted_gradient(point: np.array):
            result = grad_function(point) / coeff_function(point)
            return result

        self.__grad_function = adjusted_gradient
        self.__unknown_count = 0
        self.__unknown_values = np.empty(0)
        self.__integrator = Integrator2D(
            kernel,
            SingleLayerInclusionTerm.NOMINAL_INTEGRATION_POINTS_PER_AXIS,
            SingleLayerInclusionTerm.SINGULARITY_INTEGRATION_POINTS_PER_AXIS,
        )

    def calculate_coefficients(self, point_info: Point2DInfo, condition: BoundaryConditionType):
        assert point_info is not None
        coefficients = np.empty(0)
        for face in self.domain.get_mesh():
            is_same_point = Utils.is_the_same_point(point_info.point, face.point, SingleLayerInclusionTerm.EPS)
            if condition in [
                BoundaryConditionType.DIRICHLET,
                BoundaryConditionType.NEUMANN,
                BoundaryConditionType.INCLUSION,
            ]:
                res = 0.0
                for index in range(Constants.PLANE_SQUARE_DIM):
                    next_index = index + 1
                    if Constants.PLANE_SQUARE_DIM == next_index:
                        next_index = 0
                    normal = Utils.get_right_side_normal(face.element[index], face.element[next_index])

                    def integral_function(point: np.array):
                        result = np.dot(normal, self.__grad_function(point))
                        return result

                    res -= self.__integrator.segment_of(
                        KernelValueType.SCALAR,
                        integral_function,
                        point_info.point,
                        face.element[index],
                        face.element[next_index],
                    )
                if BoundaryConditionType.INCLUSION == condition and is_same_point:
                    res -= 1.0
                coefficients = np.append(coefficients, [res])
            elif BoundaryConditionType.ROBIN == condition:
                raise NotImplementedError("Not implemented yet")
            else:
                raise AttributeError("Not supported boundary element type")
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
        for face in self.domain.get_mesh():
            face_value = 0.0
            for index in range(Constants.PLANE_SQUARE_DIM):
                next_index = index + 1
                if Constants.PLANE_SQUARE_DIM == next_index:
                    next_index = 0
                normal = Utils.get_right_side_normal(face.element[index], face.element[next_index])

                def integral_function(point: np.array):
                    result = np.dot(normal, self.__grad_function(point))
                    return result

                face_value -= self.__integrator.segment_of(
                    KernelValueType.SCALAR,
                    integral_function,
                    point,
                    face.element[index],
                    face.element[next_index],
                )
            face_value *= self.__unknown_values[unknown_index]
            unknown_index += 1
            result += face_value
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
                if self.__domain.is_corner_point(point):
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
            if isinstance(point_info.type, list):
                conditions = point_info.type
            else:
                conditions = [point_info.type]
            for condition in conditions:
                right_side_value = 0.0
                matrix_row = np.empty(0)
                for term in self.__expression:
                    coefficients, value = term.calculate_coefficients(point_info, condition)
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

    def print_stats(self):
        print("Domains stats...")
        self.__domain.print_stats()

    def plot(self):
        print("Visualizing data...")

        x_min = min(point_info.element[0][0] for point_info in self.__domain.get_border())
        y_min = min(point_info.element[0][1] for point_info in self.__domain.get_border())
        x_max = max(point_info.element[0][0] for point_info in self.__domain.get_border())
        y_max = max(point_info.element[0][1] for point_info in self.__domain.get_border())

        x_min += GlobalSettings.CHART_SKIP_BORDER_WIDTH
        y_min += GlobalSettings.CHART_SKIP_BORDER_WIDTH
        x_max -= GlobalSettings.CHART_SKIP_BORDER_WIDTH
        y_max -= GlobalSettings.CHART_SKIP_BORDER_WIDTH

        if GlobalSettings.PLOT_DETAILS:
            functions = []
            for index in range(0, len(self.__expression)):
                functions.append(lambda point, _index=index: self.__expression[_index].value(point))
            Utils.plot([x_min, x_max], [y_min, y_max], functions, self.__domain)
        else:
            if GlobalSettings.PLOT_ERROR and self.__analytical_solution is not None:
                if GlobalSettings.PLOT_ERROR_CONTOURS:
                    Utils.plot_contour(
                        [x_min, x_max],
                        [y_min, y_max],
                        lambda point: np.abs(self.solution_value(point) - self.__analytical_solution(point)),
                        self.__domain,
                    )
                else:
                    Utils.plot(
                        [x_min, x_max],
                        [y_min, y_max],
                        [lambda point: np.abs(self.solution_value(point) - self.__analytical_solution(point))],
                        self.__domain,
                    )
                if GlobalSettings.ERROR_TO_CSV:
                    x_values = np.linspace(x_min, x_max, num=GlobalSettings.ERROR_CSV_STEPS)
                    y_values = np.linspace(y_min, y_max, num=GlobalSettings.ERROR_CSV_STEPS)

                    header = ["/"]
                    for y in y_values:
                        header.append(y)

                    data = []
                    for x in x_values:
                        line = [x]
                        for y in y_values:
                            point = np.array([x, y])
                            value = abs(self.solution_value(point) - self.__analytical_solution(point))
                            line.append(f"{value:.3}")
                        data.append(line)

                    file_path = "/home/user/workspace/project/ngsolve/solvers/samples/manual/cobem.csv"

                    with open(file_path, mode="w", newline="") as file:
                        writer = csv.writer(file)
                        writer.writerow(header)
                        writer.writerows(data)

                    print(f"Data written to {file_path}")

            else:
                if GlobalSettings.PLOT_ERROR_CONTOURS:
                    Utils.plot_contour(
                        [x_min, x_max],
                        [y_min, y_max],
                        self.solution_value,
                        self.__domain,
                    )
                else:
                    Utils.plot([x_min, x_max], [y_min, y_max], [self.solution_value], self.__domain)

        print("Done!")


class Samples:
    class Laplace:

        class BEM:

            @staticmethod
            def init_dirichlet_convex():
                print("BEM for Dirichlet problem for Laplace equation in convex shape...")

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return point[0] * point[0] - point[1] * point[1]

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = PlainDomain2D(
                    [
                        path.Path([(-0.5, -1.0), (0.5, -1.0)]),
                        path.Path([(0.5, -1.0), (0.5, 1.0)]),
                        path.Path([(0.5, 1.0), (-0.5, 1.0)]),
                        path.Path([(-0.5, 1.0), (-0.5, -1.0)]),
                    ],
                    [BoundaryConditionType.DIRICHLET] * 4,
                    [dirichlet_boundary_value] * 4,
                )

                expression = [
                    DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
                    SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
                ]

                problem = Problem(ProblemSolverType.BEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_dirichlet_single_inclusion_square():
                print("BEM for Dirichlet problem for Poisson equation with single inclusion...")

                print("Define boundary conditions...")

                INCLUSION_SIZE = 0.2
                INCLUSION_CENTER_X = 0.5
                INCLUSION_CENTER_Y = 0.5
                INCLUSION_RADIUS = INCLUSION_SIZE / 2.0 * np.sqrt(2.0)
                K_MAX = 10
                K_MIN = 1

                def boundary_value(point: np.array):
                    return 2.0 * point[1]

                inclusion_domain = SquareDomain2D(
                    np.array([INCLUSION_CENTER_X - INCLUSION_SIZE / 2.0, INCLUSION_CENTER_Y - INCLUSION_SIZE / 2.0]),
                    np.array([INCLUSION_CENTER_X + INCLUSION_SIZE / 2.0, INCLUSION_CENTER_Y + INCLUSION_SIZE / 2.0]),
                    [BoundaryConditionType.INCLUSION] * 4,
                    [Utils.constant_one()] * 4,
                    GlobalSettings.INCLUSION_ELEMENTS_COUNT,
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
                    result = K_MIN
                    distance_from_center = (
                        (point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2
                    ) / INCLUSION_RADIUS**2
                    if distance_from_center < 1.0:
                        result += K_MAX * (1 - distance_from_center)

                    return result

                def thermal_conductivity_gradient(point: np.array):
                    result = np.array([0.0, 0.0])
                    if inclusion_domain.is_point_inside_domain(point):
                        result[0] = K_MAX * 2.0 * (INCLUSION_CENTER_X - point[0]) / INCLUSION_RADIUS**2
                        result[1] = K_MAX * 2.0 * (INCLUSION_CENTER_Y - point[1]) / INCLUSION_RADIUS**2
                    return result

                expression = [
                    SingleLayerInclusionTerm(
                        Laplace2DKernel(),
                        inclusion_domain,
                        thermal_conductivity,
                        thermal_conductivity_gradient,
                        1.0,
                    ),
                    DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
                    SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
                ]

                problem = Problem(ProblemSolverType.BEM, expression, domain)
                return problem

            @staticmethod
            def init_dirichlet_single_inclusion_circle_in_square():
                print("BEM for Dirichlet problem for Laplace equation with single circular inclusion...")

                print("Define boundary conditions...")

                INCLUSION_CENTER = np.array([0.5, 0.5])
                INCLUSION_RADIUS = 0.2
                K_MAX = 10
                K_MIN = 1

                def boundary_value(point: np.array):
                    return 2.0 * point[1]

                inclusion_domain = GmshDomain2D(
                    MeshLoader.generate_mesh_from_file(
                        MeshLoader.UNIT_CIRCLE_FILE, GlobalSettings.INCLUSION_ELEMENTS_MAX_SIZE, INCLUSION_RADIUS
                    ),
                    MeshLoader.order_conditions(
                        {
                            MeshLoader.LEFT: (BoundaryConditionType.INCLUSION, Utils.constant_one()),
                            MeshLoader.RIGHT: (BoundaryConditionType.INCLUSION, Utils.constant_one()),
                        }
                    ),
                    INCLUSION_CENTER,
                )

                domain = GmshDomain2D(
                    MeshLoader.generate_mesh_from_file(
                        MeshLoader.UNIT_SQUARE_FILE, GlobalSettings.BORDER_ELEMENT_MAX_SIZE
                    ),
                    MeshLoader.order_conditions(
                        {
                            MeshLoader.TOP: (BoundaryConditionType.DIRICHLET, boundary_value),
                            MeshLoader.LEFT: (BoundaryConditionType.DIRICHLET, boundary_value),
                            MeshLoader.RIGHT: (BoundaryConditionType.DIRICHLET, boundary_value),
                            MeshLoader.BOTTOM: (BoundaryConditionType.DIRICHLET, boundary_value),
                        }
                    ),
                    np.array([0.5, 0.5]),
                    [inclusion_domain],
                )

                def thermal_conductivity(point: np.array):
                    result = K_MIN
                    distance_from_center = Utils.squared_distance(point, INCLUSION_CENTER) / INCLUSION_RADIUS**2
                    if distance_from_center < 1.0:
                        result += K_MAX * (1 - distance_from_center)
                    return result

                def thermal_conductivity_gradient(point: np.array):
                    result = np.array([0.0, 0.0])
                    if Utils.distance(point, INCLUSION_CENTER) < INCLUSION_RADIUS:
                        result = K_MAX * 2.0 * (INCLUSION_CENTER - point) / INCLUSION_RADIUS**2
                    return result

                expression = [
                    SingleLayerInclusionTerm(
                        Laplace2DKernel(),
                        inclusion_domain,
                        thermal_conductivity,
                        thermal_conductivity_gradient,
                        1.0,
                    ),
                    DoubleLayerBoundaryTerm(Laplace2DKernel(), domain),
                    SingleLayerBoundaryTerm(Laplace2DKernel(), domain, -1),
                ]

                problem = Problem(ProblemSolverType.BEM, expression, domain)
                return problem

        class CoBEM:

            @staticmethod
            def init_dirichlet_hexagon():
                print("CoBEM for Dirichlet problem for Laplace equation in hexagon...")

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return (point[0] + 2) * (point[0] + 2) - (point[1] + 1) * (point[1] + 1)

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = PlainDomain2D(
                    [
                        path.Path([(-2.0, 0.0), (-1.0, -1.0), (1.0, -1.0), (2.0, 0.0), (1.0, 1.0), (-1.0, 1.0)]),
                        path.Path([(-1.0, 1.0), (-2.0, 0.0)]),
                    ],
                    [BoundaryConditionType.DIRICHLET] * 2,
                    [dirichlet_boundary_value] * 2,
                )

                expression = [
                    SingleLayerCoBEMTerm(Laplace2DKernel(), domain, 1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_dirichlet_convex():
                print("CoBEM for Dirichlet problem for Laplace equation in convex shape...")

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return point[0] * point[0] - point[1] * point[1]

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = PlainDomain2D(
                    [
                        path.Path([(-0.5, -1.0), (0.5, -1.0)]),
                        path.Path([(0.5, -1.0), (0.5, 1.0)]),
                        path.Path([(0.5, 1.0), (-0.5, 1.0)]),
                        path.Path([(-0.5, 1.0), (-0.5, -1.0)]),
                    ],
                    [BoundaryConditionType.DIRICHLET] * 4,
                    [dirichlet_boundary_value] * 4,
                )

                expression = [
                    SingleLayerCoBEMTerm(Laplace2DKernel(), domain, 1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

    class Poisson:

        class BEM:

            @staticmethod
            def init_dirichlet_square():
                print("BEM for Dirichlet problem for Poisson equation...")

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return 2 * (point[0] - 0.5) ** 2 + 2 * (point[1] - 0.5) ** 2

                def boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = GmshDomain2D(
                    MeshLoader.generate_mesh_from_file(
                        MeshLoader.UNIT_SQUARE_FILE, GlobalSettings.BORDER_ELEMENT_MAX_SIZE
                    ),
                    MeshLoader.order_conditions(
                        {
                            MeshLoader.TOP: (BoundaryConditionType.DIRICHLET, boundary_value),
                            MeshLoader.LEFT: (BoundaryConditionType.DIRICHLET, boundary_value),
                            MeshLoader.RIGHT: (BoundaryConditionType.DIRICHLET, boundary_value),
                            MeshLoader.BOTTOM: (BoundaryConditionType.DIRICHLET, boundary_value),
                        }
                    ),
                    np.array([0.5, 0.5]),
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

            @staticmethod
            def init_neumann_square():
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

            @staticmethod
            def init_robin_square():
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

            @staticmethod
            def init_dirichlet_convex():
                print("BEM for Dirichlet problem for Poisson equation in convex...")

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return 2 * (point[0] - 0.5) ** 2 + 2 * (point[1] - 0.5) ** 2

                def boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = PlainDomain2D(
                    [
                        path.Path([(-0.5, -1.0), (0.5, -1.0)]),
                        path.Path([(0.5, -1.0), (0.5, 1.0)]),
                        path.Path([(0.5, 1.0), (-0.5, 1.0)]),
                        path.Path([(-0.5, 1.0), (-0.5, -1.0)]),
                    ],
                    [BoundaryConditionType.DIRICHLET] * 4,
                    [boundary_value] * 4,
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

        class CoBEM:

            @staticmethod
            def init_dirichlet_square():
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
                    SingleLayerCoBEMTerm(Laplace2DKernel(), domain),
                    SingleLayerVolumeCoBEMTerm(Laplace2DKernel(), domain, heat_source_function, -1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_neumann_square():
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
                    SingleLayerCoBEMTerm(Laplace2DKernel(), domain),
                    SingleLayerVolumeCoBEMTerm(Laplace2DKernel(), domain, heat_source_function, -1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_robin_square():
                print("CoBEM for Robin problem for Poisson equation...")

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
                    SingleLayerCoBEMTerm(Laplace2DKernel(), domain),
                    SingleLayerVolumeCoBEMTerm(Laplace2DKernel(), domain, heat_source_function, -1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_dirichlet_single_inclusion_square():
                print("CoBEM for Dirichlet problem for Poisson equation with single inclusion...")

                print("Define boundary conditions...")

                INCLUSION_SIZE = 0.2
                INCLUSION_CENTER_X = 0.5
                INCLUSION_CENTER_Y = 0.5
                INCLUSION_RADIUS = INCLUSION_SIZE / 2.0 * np.sqrt(2.0)
                K_MAX = 10
                K_MIN = 1

                def boundary_value(point: np.array):
                    return 2.0 * point[1]

                inclusion_domain = SquareDomain2D(
                    np.array([INCLUSION_CENTER_X - INCLUSION_SIZE / 2.0, INCLUSION_CENTER_Y - INCLUSION_SIZE / 2.0]),
                    np.array([INCLUSION_CENTER_X + INCLUSION_SIZE / 2.0, INCLUSION_CENTER_Y + INCLUSION_SIZE / 2.0]),
                    [BoundaryConditionType.INCLUSION] * 4,
                    [Utils.constant_one()] * 4,
                    GlobalSettings.INCLUSION_ELEMENTS_COUNT,
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
                    result = K_MIN
                    distance_from_center = (
                        (point[0] - INCLUSION_CENTER_X) ** 2 + (point[1] - INCLUSION_CENTER_Y) ** 2
                    ) / INCLUSION_RADIUS**2
                    if distance_from_center < 1.0:
                        result += K_MAX * (1 - distance_from_center)

                    return result

                def thermal_conductivity_gradient(point: np.array):
                    result = np.array([0.0, 0.0])
                    if inclusion_domain.is_point_inside_domain(point):
                        result[0] = K_MAX * 2.0 * (INCLUSION_CENTER_X - point[0]) / INCLUSION_RADIUS**2
                        result[1] = K_MAX * 2.0 * (INCLUSION_CENTER_Y - point[1]) / INCLUSION_RADIUS**2
                    return result

                print("Define heat source function...")

                expression = [
                    SingleLayerInclusionTerm(
                        Laplace2DKernel(),
                        inclusion_domain,
                        thermal_conductivity,
                        thermal_conductivity_gradient,
                        1.0,
                    ),
                    SingleLayerCoBEMTerm(Laplace2DKernel(), domain),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain)
                return problem

    class Pennes:

        class BEM:

            # Article https://pubmed.ncbi.nlm.nih.gov/1522731/
            @staticmethod
            def init_dirichlet_square():
                print("BEM for Dirichlet problem for Pennes equation...")

                K_SQUARE_CONSTANT = 2.0

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return np.sinh(point[0] + point[1])

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = SquareDomain2D(
                    np.array([0.0, 0.0]),
                    np.array([1.0, 1.0]),
                    [BoundaryConditionType.DIRICHLET] * 4,
                    [dirichlet_boundary_value] * 4,
                    GlobalSettings.BORDER_ELEMENTS_COUNT,
                )

                expression = [
                    DoubleLayerBoundaryTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain),
                    SingleLayerBoundaryTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain, -1),
                ]

                problem = Problem(ProblemSolverType.BEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_neumann_square():
                print("BEM for Neumann problem for Pennes equation...")

                K_SQUARE_CONSTANT = 1.0

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return np.exp(point[1])

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                def neumann_boundary_right_value(point: np.array):
                    return 0.0

                def neumann_boundary_left_value(point: np.array):
                    return 0.0

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

                expression = [
                    DoubleLayerBoundaryTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain),
                    SingleLayerBoundaryTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain, -1),
                ]

                problem = Problem(ProblemSolverType.BEM, expression, domain, analytical_solution)
                return problem

        class CoBEM:

            @staticmethod
            def init_dirichlet_square():
                print("CoBEM for Dirichlet problem for Pennes equation...")

                K_SQUARE_CONSTANT = 2.0

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return np.sinh(point[0] + point[1])

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = SquareDomain2D(
                    np.array([0.0, 0.0]),
                    np.array([1.0, 1.0]),
                    [BoundaryConditionType.DIRICHLET] * 4,
                    [dirichlet_boundary_value] * 4,
                    GlobalSettings.BORDER_ELEMENTS_COUNT,
                )

                expression = [
                    SingleLayerCoBEMTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain, 1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_neumann_square():
                print("CoBEM for Neumann problem for Pennes equation...")

                K_SQUARE_CONSTANT = 1.0

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return np.exp(point[1])

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                def neumann_boundary_right_value(point: np.array):
                    return 0.0

                def neumann_boundary_left_value(point: np.array):
                    return 0.0

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

                expression = [
                    SingleLayerCoBEMTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain, 1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem

            @staticmethod
            def init_dirichlet_hexagon():
                print("CoBEM for Dirichlet problem for Pennes equation in hexagon...")

                K_SQUARE_CONSTANT = 0.25

                print("Define boundary conditions...")

                def analytical_solution(point: np.array):
                    return np.sinh(0.5 * (point[1] + 1))

                def dirichlet_boundary_value(point: np.array):
                    return analytical_solution(point)

                domain = PlainDomain2D(
                    [
                        path.Path([(-3.0, 0.0), (-2.0, -1.0), (2.0, -1.0), (3.0, 0.0), (2.0, 1.0), (-2.0, 1.0)]),
                        path.Path([(-2.0, 1.0), (-3.0, 0.0)]),
                    ],
                    [BoundaryConditionType.DIRICHLET] * 2,
                    [dirichlet_boundary_value] * 2,
                )

                expression = [
                    SingleLayerCoBEMTerm(Pennes2DKernel(K_SQUARE_CONSTANT), domain, 1),
                ]

                problem = Problem(ProblemSolverType.COBEM, expression, domain, analytical_solution)
                return problem


if "__main__" == __name__:
    problem = None
    if 1 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.BEM.init_dirichlet_square()
    elif 2 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.BEM.init_neumann_square()
    elif 3 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.BEM.init_robin_square()
    elif 4 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Laplace.BEM.init_dirichlet_single_inclusion_square()
    elif 5 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.CoBEM.init_dirichlet_square()
    elif 6 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.CoBEM.init_neumann_square()
    elif 7 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.CoBEM.init_robin_square()
    elif 8 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.CoBEM.init_dirichlet_single_inclusion_square()
    elif 9 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Pennes.BEM.init_dirichlet_square()
    elif 10 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Pennes.BEM.init_neumann_square()
    elif 11 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Pennes.CoBEM.init_dirichlet_square()
    elif 12 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Pennes.CoBEM.init_neumann_square()
    elif 13 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Pennes.CoBEM.init_dirichlet_hexagon()
    elif 14 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Laplace.CoBEM.init_dirichlet_hexagon()
    elif 15 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Laplace.BEM.init_dirichlet_convex()
    elif 16 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Laplace.CoBEM.init_dirichlet_convex()
    elif 17 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Poisson.BEM.init_dirichlet_convex()
    elif 18 == GlobalSettings.EXAMPLE_TYPE:
        pass
        # problem = Samples.Poisson.CoBEM.init_dirichlet_convex()
    elif 19 == GlobalSettings.EXAMPLE_TYPE:
        problem = Samples.Laplace.BEM.init_dirichlet_single_inclusion_circle_in_square()

    assert problem is not None

    problem.calculate()
    problem.print_stats()
    problem.plot()
