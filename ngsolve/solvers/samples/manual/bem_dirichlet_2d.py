#!/usr/bin/env python3

import numpy as np

class ExpressionTerm:
    pass

class BoundaryCondition:
    def __init__(self):
        pass
    
    def value(self, point: np.array):
        raise NotImplementedError("Call to abstract method")

class BottomDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__()
        self.__border_elements = np.array([[1.0, 0.0], [0.0, 1.0]])

    def value(self, point: np.array):
        return 0.0

class LeftDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__()
        self.__border_elements = np.array([[0.0, 0.0], [1.0, 0.0]])

    def value(self, point: np.array):
        return 2.0 * point[1]

class RightDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__()
        self.__border_elements = np.array([[1.0, 0.0], [1.0, 1,0]])

    def value(self, point: np.array):
        return 2.0 * point[1]

class TopDirichletCondition(BoundaryCondition):
    def __init__(self):
        super().__init__()
        self.__border_elements = np.array([[1.0, 0.0], [1.0, 1.0]])

    def value(self, point: np.array):
        return 2.0

class ExpressionTerm:
    pass

class SingleLayerBoundaryTerm(ExpressionTerm):
    pass

class DoubleLayerBoundaryTerm(ExpressionTerm):
    pass

class SingleLayerVolumeTerm(ExpressionTerm):
    pass

print("Define boundary conditions...")

boundary_conditions = [BottomDirichletCondition(), LeftDirichletCondition(), RightDirichletCondition(), 
                       TopDirichletCondition()]

expression = [SingleLayerBoundaryTerm(), DoubleLayerBoundaryTerm(), SingleLayerVolumeTerm()]

unknown_count = 10

print("Create SLA...")

matrix = None
right_side = None

for boundary in boundary_conditions:
    for point in boundary.get_boundary_points():
        matrix_row = np.zeros(unknown_count)
        right_side_value = 0.0
        for term in expression:
            term.calculate_coefficients(matrix_row)
            right_side += term.calculate_right_side()
        if matrix is None:
            matrix = matrix_row
        else:
            matrix = np.vstack((matrix, matrix_row))
        if right_side is None:
            right_side = np.array([right_side_value])
        else:
            right_side = np.append(right_side, [right_side_value])
        

print("Solving SLA...")

print("Setting data back...")

print("Visualizing data...")

print("Done!")
