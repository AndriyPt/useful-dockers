from netgen.occ import *
from ngsolve import *
from ngsolve import VTKOutput
from typing import List

E_TISSUE = 1060 # tissue
NU_TISSUE = 0.31

MU_TISSUE  = E_TISSUE / 2 / (1 + NU_TISSUE)
LAMBDA_TISSUE = E_TISSUE * NU_TISSUE / ((1 + NU_TISSUE) * (1 - 2 * NU_TISSUE))

E_TUMOR = 5460 + 3180 # tumor
NU_TUMOR = 0.26

MU_TUMOR  = E_TUMOR / 2 / (1 + NU_TUMOR)
LAMBDA_TUMOR = E_TUMOR * NU_TUMOR / ((1 + NU_TUMOR) * (1 - 2 * NU_TUMOR))

DIMENSIONS = 2

BORDER_TOP = "gamma_top"
BORDER_PRESS = "gamma_pressure"
BORDER_SIDE = "gamma_side"
BORDER_BOTTOM = "gamma_bottom"

DOMAIN_TISSUE = "tissue"
DOMAIN_TUMOR = "tumor"

PARAM_A = 0.04 # height (m)
PARAM_B = 0.16 # CO2 area height (m)
PARAM_C = 0.4 # width (m)

TOOL_RADIUS = 0.0025 # m
TOOL_DISPLACEMENT = -0.001 # m

LARGE_MAXH=0.0025

class TumorParams:
    def __init__(self, x:float = 0.0, y:float = 0.0, radius:float = 0.0, column_name:str = "y"):
        self.col_name = column_name
        self.radius = radius
        self.center_x = x
        self.center_y = y

def calculate(params: TumorParams, x_values:List[float], y_value:float, visualize:bool = False):

    # generate a triangular mesh
    left_rect = Rectangle(PARAM_C / 2.0 - TOOL_RADIUS, PARAM_A).Face()
    left_rect.edges.Min(X).name = BORDER_SIDE
    left_rect.edges.Max(Y).name = BORDER_TOP
    left_rect.edges.Min(Y).name = BORDER_BOTTOM
    left_rect.faces.name = DOMAIN_TISSUE

    middle_rect = MoveTo(PARAM_C / 2.0 - TOOL_RADIUS, 0).Rectangle(2.0 * TOOL_RADIUS, PARAM_A).Face()
    middle_rect.edges.Max(Y).name = BORDER_PRESS
    middle_rect.edges.Min(Y).name = BORDER_BOTTOM
    middle_rect.faces.name = DOMAIN_TISSUE

    right_rect = MoveTo(PARAM_C / 2.0 + TOOL_RADIUS, 0).Rectangle(PARAM_C / 2.0 - TOOL_RADIUS, PARAM_A).Face()
    right_rect.edges.Max(X).name = BORDER_SIDE
    right_rect.edges.Max(Y).name = BORDER_TOP
    right_rect.edges.Min(Y).name = BORDER_BOTTOM
    right_rect.faces.name = DOMAIN_TISSUE

    whole_body = Glue([left_rect, middle_rect, right_rect])
    whole_body.faces.name = DOMAIN_TISSUE

    tumor_cross_section = Circle((params.center_x, params.center_y), params.radius).Face()
    tumor_cross_section.faces.name = DOMAIN_TUMOR

    tissue_shape = whole_body - tumor_cross_section 
    tissue_shape.faces.name = DOMAIN_TISSUE

    shape = Glue([tissue_shape, tumor_cross_section])

    geo = OCCGeometry(shape, dim = DIMENSIONS)

    mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

    print("Boundaries: ", mesh.GetBoundaries())
    print("Materials: ", mesh.GetMaterials())

    if visualize:
        Draw(mesh)

    lambda_coef = mesh.MaterialCF({
        DOMAIN_TISSUE: LAMBDA_TISSUE,
        DOMAIN_TUMOR: (LAMBDA_TUMOR - LAMBDA_TISSUE) * cos(0.5 * pi / params.radius ** 2 * ((x - params.center_x) ** 2 + 
                                                    (y - params.center_y) ** 2)) + LAMBDA_TISSUE,
        }, 
        default = 0)

    mu_coef = mesh.MaterialCF({
        DOMAIN_TISSUE: LAMBDA_TISSUE,
        DOMAIN_TUMOR: (LAMBDA_TUMOR - LAMBDA_TISSUE) * cos(0.5 * pi / params.radius ** 2 * ((x - params.center_x) ** 2 + 
                                                    (y - params.center_y) ** 2)) + LAMBDA_TISSUE,
        }, 
        default = 0)


    def Stress(strain):
        return 2 * mu_coef * strain + lambda_coef * Trace(strain) * Id(DIMENSIONS)   

    fes = VectorH1(mesh, order=3, dirichlet=BORDER_BOTTOM + "|" + BORDER_SIDE + "|" + BORDER_PRESS)

    # Dirichlet conditions
    dirichlet_conditions = mesh.BoundaryCF({BORDER_PRESS: (0, TOOL_DISPLACEMENT)}, default = (0,0))
    dirichlet_gfu = GridFunction(fes)
    dirichlet_gfu.Set(dirichlet_conditions, BND)

    u,v = fes.TnT()
    gfu = GridFunction(fes)

    with TaskManager():
        a = BilinearForm(InnerProduct(Stress(Sym(Grad(u))), Sym(Grad(v))).Compile()*dx)
        pre = Preconditioner(a, "bddc")
        a.Assemble()

    from ngsolve.krylovspace import CGSolver
    inv = CGSolver(a.mat, pre, tol=1e-8)

    # Approach for nonhomogeneous Dirichlet boundary condition
    # https://docu.ngsolve.org/release/i-tutorials/unit-1.3-dirichlet/dirichlet.html
    f = LinearForm(fes) # Used to create vector of needed dimension 
    r = f.vec.CreateVector()
    r.data = inv * dirichlet_gfu.vec
    gfu.vec.data = dirichlet_gfu.vec.data - inv * r

    with TaskManager():
        fesstress = MatrixValued(H1(mesh, order=3), symmetric=True)
        gfstress = GridFunction(fesstress)
        gfstress.Interpolate(Stress(Sym(Grad(gfu))))


    normal = CoefficientFunction((0.0, 1.0))

    traction = gfstress * normal 

    print ("Force applied (Newtons):", TOOL_RADIUS * Integrate(traction, mesh, definedon=mesh.Boundaries(BORDER_PRESS)))
    
    if visualize:
        Draw(gfu, mesh, "Deformations")
        Draw(traction, mesh, "Traction")

    if visualize:
        vtk = VTKOutput(ma=mesh,
                        coefs=[gfu, traction],
                        names = ["displacement", "traction"],
                        filename="/tmp/elasticity_2d_result",
                        subdivision=3)
        vtk.Do()
    
    result = []

    for x_point in x_values:
        point = mesh(x_point, y_value) 
        value_traction = traction(point)[1] # Get only the second component of the value
        result.append(value_traction)

    return result
    

def generate_data(tumor_params: List[TumorParams], debug_index:int = -1, visualize:bool = False, 
                  write_to_csv:bool = True):

    count = 100
    step = 4.0 * 2.0 * TOOL_RADIUS / count
    x_values = [PARAM_C / 2.0 - 4.0 * TOOL_RADIUS + x * step for x in range(count + 1)]
    y_value = PARAM_A

    data = []

    if (debug_index < 0):
        for param in tumor_params:
            result = calculate(param, x_values, y_value, visualize)
            data.append(result)
    else:
        calculate(tumor_params[debug_index], x_values, y_value, visualize)

    if write_to_csv and debug_index < 0:
        with open("/home/user/workspace/project/docs/dissertation/data/elasticity_2d_border_data.csv", "w") as data_csv:
            header = "x"
            for param in tumor_params:
                header += ",{}".format(param.col_name)
            data_csv.write("{}\n".format(header))

            for i in range(len(x_values)):
                line = "{}".format(x_values[i])
                for j in range(len(tumor_params)):
                    line += ",{}".format(data[j][i])
                data_csv.write("{}\n".format(line))        

def main():
    tumor_radius = 0.0025 # m
    larger_tumor_radius = 0.005 # m
    extra_depth = 0.001 # m

    tumor_params = [
        TumorParams(PARAM_C / 2.0, PARAM_A - tumor_radius, tumor_radius, "y_d_0025_c_0"),
        TumorParams(PARAM_C / 2.0, PARAM_A - tumor_radius - extra_depth, tumor_radius, "y_d_0025_c_001"),
        TumorParams(PARAM_C / 2.0, PARAM_A - larger_tumor_radius, larger_tumor_radius, "y_d_005_c_0"),
        TumorParams(PARAM_C / 2.0, PARAM_A - larger_tumor_radius - extra_depth, larger_tumor_radius, "y_d_005_c_001"),
        ]
    
    generate_data(tumor_params, debug_index=-1, visualize=False, write_to_csv=True)
    # generate_data(tumor_params, debug_index=0, visualize=True, write_to_csv=False)

    print("Done!")

main()

