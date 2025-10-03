from bem_dirichlet_2d import Kernel, Integrator2D, HexagonalDomain2D, BoundaryConditionType, GlobalSettings
import numpy as np

FLOAT_VALUE_EPSILON = 0.001


def assert_floats_are_equal(expected: float, actual: float):
    assert abs(expected - actual) < FLOAT_VALUE_EPSILON, "Expected {} but received {}".format(expected, actual)


class IdentityKernel(Kernel):
    def value(self, point_x: np.array, point_y: np.array):
        return 1.0

    def grad(self, point_x: np.array, point_y: np.array):
        return np.zeros(2)


def identity_value(point: np.array):
    return 1.0


def test_segment_integration():
    NOMINAL_POINTS = 4
    SINGULARITY_POINTS = 6
    integrator = Integrator2D(IdentityKernel(), NOMINAL_POINTS, SINGULARITY_POINTS)

    result = integrator.segment(identity_value, np.array([0.5, 0.5]), np.array([0.0, 0.0]), np.array([1.0, 0.0]))
    assert_floats_are_equal(1.0, result)

    result = integrator.segment(identity_value, np.array([0.5, 0.5]), np.array([0.0, 0.0]), np.array([2.0, 0.0]))
    assert_floats_are_equal(2.0, result)

    result = integrator.segment(identity_value, np.array([0.5, 0.5]), np.array([0.0, 0.0]), np.array([5.0, 0.0]))
    assert_floats_are_equal(5.0, result)


def test_square_integration():
    NOMINAL_POINTS = 4
    SINGULARITY_POINTS = 6
    integrator = Integrator2D(IdentityKernel(), NOMINAL_POINTS, SINGULARITY_POINTS)

    unit_square = np.array([np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([1.0, 1.0]), np.array([0.0, 1.0])])
    result = integrator.square(
        identity_value,
        np.array([0.5, 0.5]),
        unit_square,
    )
    assert_floats_are_equal(1.0, result)

    result = integrator.square(identity_value, np.array([0.5, 0.5]), 2.0 * unit_square)
    assert_floats_are_equal(4.0, result)

    result = integrator.square(identity_value, np.array([0.5, 0.5]), 5.0 * unit_square)
    assert_floats_are_equal(25.0, result)


def test_trapezoid_integration():
    NOMINAL_POINTS = 4
    SINGULARITY_POINTS = 6
    integrator = Integrator2D(IdentityKernel(), NOMINAL_POINTS, SINGULARITY_POINTS)

    result = integrator.trapezoid(
        identity_value,
        np.array([0.5, 0.5]),
        np.array([np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([1.0, 2.0]), np.array([0.0, 1.0])]),
    )
    assert_floats_are_equal(1.5, result)

    result = integrator.trapezoid(
        identity_value,
        np.array([0.5, 0.5]),
        np.array([np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([1.0, 1.0]), np.array([0.0, 2.0])]),
    )
    assert_floats_are_equal(1.5, result)

    result = integrator.trapezoid(
        identity_value,
        np.array([0.5, 0.5]),
        np.array([np.array([0.0, 0.0]), np.array([2.0, 0.0]), np.array([2.0, 2.0]), np.array([0.0, 4.0])]),
    )
    assert_floats_are_equal(6, result)

    result = integrator.trapezoid(
        identity_value,
        np.array([0.5, 0.5]),
        np.array([np.array([-3.0, 0.0]), np.array([3.0, 0.0]), np.array([3.0, 3.0]), np.array([0.0, 3.0])]),
    )
    assert_floats_are_equal(13.5, result)

    result = integrator.trapezoid(
        identity_value,
        np.array([0.5, 0.5]),
        np.array([np.array([0.0, 0.0]), np.array([3.0, 0.0]), np.array([2.0, 2.0]), np.array([0.0, 2.0])]),
    )
    assert_floats_are_equal(5.0, result)

    result = integrator.trapezoid(
        identity_value,
        np.array([0.5, 0.5]),
        np.array([np.array([0.0, 0.0]), np.array([3.0, 0.0]), np.array([2.0, 2.0]), np.array([1.0, 2.0])]),
    )
    assert_floats_are_equal(4.0, result)


def test_hexagonal_domain():

    GlobalSettings.COBORDER_DEPTH = 1.0

    domain = HexagonalDomain2D(
        np.array([-1.0, 0.0]),
        np.array([1.0, 2.0]),
        4.0,
        [BoundaryConditionType.DIRICHLET] * 6,
        [lambda point: 0.0] * 6,
        1
    )

    assert 6 == len(domain.get_edge_points())

    assert_floats_are_equal(-1.0, domain.get_edge_points()[0][0])
    assert_floats_are_equal(0.0, domain.get_edge_points()[0][1])

    assert_floats_are_equal(1.0, domain.get_edge_points()[1][0])
    assert_floats_are_equal(0.0, domain.get_edge_points()[1][1])

    assert_floats_are_equal(2.0, domain.get_edge_points()[2][0])
    assert_floats_are_equal(1.0, domain.get_edge_points()[2][1])

    assert_floats_are_equal(1.0, domain.get_edge_points()[3][0])
    assert_floats_are_equal(2.0, domain.get_edge_points()[3][1])

    assert_floats_are_equal(-1.0, domain.get_edge_points()[4][0])
    assert_floats_are_equal(2.0, domain.get_edge_points()[4][1])

    assert_floats_are_equal(-2.0, domain.get_edge_points()[5][0])
    assert_floats_are_equal(1.0, domain.get_edge_points()[5][1])

    assert 6 == len(domain.get_border())

    assert_floats_are_equal(0.0, domain.get_border()[0].point[0])
    assert_floats_are_equal(0.0, domain.get_border()[0].point[1])

    assert_floats_are_equal(1.5, domain.get_border()[1].point[0])
    assert_floats_are_equal(0.5, domain.get_border()[1].point[1])

    assert 6 == len(domain.get_coborder())

def main():
    test_segment_integration()
    test_square_integration()
    test_trapezoid_integration()
    test_hexagonal_domain()
    print("All tests have passed!")


if "__main__" == __name__:
    main()
