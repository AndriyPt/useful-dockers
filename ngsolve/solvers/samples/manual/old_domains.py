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
        self._init_coborder_polygon(0.5 * (bottom_left_point + top_right_point))

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

            item.element = Utils.order_quad_ccw(elements)
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

    def print_stats(self):
        print("---")
        print("SquareDomain stats:")
        print("---")
        print(f"Border elements count: {len(self.__border)}")
        print(f"Squares count: {len(self.__mesh)}")
        print("---")
        for domain in self.get_subdomains():
            domain.print_stats()


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

class Samples:
    class Laplace:

        class BEM:

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

        class CoBEM:

            @staticmethod
            def init_dirichlet_single_inclusion_square():
                print("CoBEM for Dirichlet problem for Laplace equation with single inclusion...")

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

