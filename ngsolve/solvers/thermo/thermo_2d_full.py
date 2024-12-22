from ngsolve import *
from netgen.occ import *
from ngsolve import VTKOutput
from typing import List

# Bladder with tumor

ROBIN_H = 8.368 # W/m^2/^K
ROBIN_T_e = 21 # ^C m, Typical OR temperature

K_TISSUE = 0.19 # W/m/^C
K_MAX_TUMOR = 0.495 # W/m/^C
K_CO2 = 0.0176 # W/m/^C

THERMAL_Q_m_TISSUE = 1.2 * 1086 # W/m^3, based on data Urinary Bladder Wall 1.21 W/kg, density 1086 kg/m**3
THERMAL_Q_m_TUMOR = THERMAL_Q_m_TISSUE * 10 # W/m^3
THERMAL_T_a = 36.8 # ^C  
DENSITY_B = 1050 # kg/m^3
THERMAL_C_b = 3617 # J/kg/^C
THERMAL_omega_TISSUE = 0.0002715 # m^3 / s / m^3
THERMAL_omega_TUMOR = 4 * THERMAL_omega_TISSUE # m^3 / s / m^3

THERMAL_W_TISSUE = DENSITY_B * THERMAL_C_b * THERMAL_omega_TISSUE
THERMAL_W_TUMOR = DENSITY_B * THERMAL_C_b * THERMAL_omega_TUMOR

BORDER_GAMMA1 = "gamma1"
BORDER_GAMMA2 = "gamma2"

DOMAIN_TISSUE = "tissue"
DOMAIN_TUMOR = "tumor"
DOMAIN_CO2 = "co2"

PARAM_A = 0.04 # height (m)
PARAM_B = 0.16 # CO2 area height (m)
PARAM_C = 0.4 # width (m)

LARGE_MAXH=0.005

class TumorParams:
    def __init__(self, x:float = 0.0, y:float = 0.0, radius:float = 0.0, column_name:str = "y"):
        self.col_name = column_name
        self.radius = radius
        self.center_x = x
        self.center_y = y


def calculate(params: TumorParams, x_values:List[float], y_value:float, visualize:bool = False):

    # generate a triangular mesh
    whole_area = Rectangle(PARAM_C, PARAM_A + PARAM_B).Face()
    whole_area.edges.name = BORDER_GAMMA2
    whole_area.edges.Min(Y).name = BORDER_GAMMA1

    bottom_area = Rectangle(PARAM_C, PARAM_A).Face()
    co2_cross_section = whole_area - bottom_area 
    co2_cross_section.faces.name = DOMAIN_CO2

    tumor_cross_section = Circle((params.center_x, params.center_y), params.radius).Face()
    tumor_cross_section.faces.name = DOMAIN_TUMOR

    tissue_cross_section = whole_area - co2_cross_section - tumor_cross_section 
    tissue_cross_section.faces.name = DOMAIN_TISSUE

    shape = Glue([tissue_cross_section, tumor_cross_section, co2_cross_section])

    geo = OCCGeometry(shape, dim = 2)
    mesh = Mesh(geo.GenerateMesh(maxh=LARGE_MAXH)).Curve(3)

    if visualize:
        print("Boundaries: ", mesh.GetBoundaries())
        print("Materials: ", mesh.GetMaterials())

        Draw(mesh)

    thermal_conductivity = mesh.MaterialCF({
        DOMAIN_TISSUE: K_TISSUE,
        DOMAIN_TUMOR: (K_MAX_TUMOR - K_TISSUE) * cos(0.5 * pi / params.radius ** 2 * (
            (x - params.center_x) ** 2 + (y - params.center_y) ** 2)) + K_TISSUE,
        DOMAIN_CO2: K_CO2,
        }, 
        default = 0)

    if visualize:
        Draw(thermal_conductivity, mesh, "Thermal Conductivity")

    heat_source = mesh.MaterialCF({ 
        DOMAIN_TISSUE: THERMAL_W_TISSUE * THERMAL_T_a + THERMAL_Q_m_TISSUE,
        DOMAIN_TUMOR: THERMAL_W_TUMOR * THERMAL_T_a + THERMAL_Q_m_TUMOR,
        },
        default = 0)

    if visualize:
        Draw(heat_source, mesh, "Heat Source")

    u_coeficient = mesh.MaterialCF({ 
        DOMAIN_TISSUE: THERMAL_W_TISSUE,
        DOMAIN_TUMOR: THERMAL_W_TUMOR,
        },
        default = 0)

    if visualize:
        Draw(u_coeficient, mesh, "u Coefficient")

    # H1-conforming finite element space
    fes = H1(mesh, order=3, dirichlet=BORDER_GAMMA1)

    # Dirichlet conditions
    dirichlet_conditions = mesh.BoundaryCF({BORDER_GAMMA1: THERMAL_T_a}, default = 0)
    dirichlet_gfu = GridFunction(fes)
    dirichlet_gfu.Set(dirichlet_conditions, BND)

    if visualize:
        Draw(dirichlet_gfu, mesh, "Partial Solution")

    # define trial- and test-functions
    u = fes.TrialFunction()
    v = fes.TestFunction()

    # the right hand side
    f = LinearForm(fes)
    f += heat_source * v * dx + ROBIN_H * ROBIN_T_e * v * ds(BORDER_GAMMA2)

    # the bilinear-form 
    a = BilinearForm(fes, symmetric=True)
    a += (thermal_conductivity * grad(u) * grad(v) + u_coeficient * u * v) * dx + ROBIN_H * u * v * ds(BORDER_GAMMA2) 

    a.Assemble()
    f.Assemble()

    # the solution field 
    gfu = GridFunction(fes)

    # Approach for nonhomogeneous Dirichlet boundary condition
    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * dirichlet_gfu.vec    

    # the solution field 
    gfu.vec.data = dirichlet_gfu.vec.data + a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r

    if visualize:
        Draw(gfu)

    flux = -grad(gfu)
    if visualize:
        Draw(flux, mesh, "Flux")

    if visualize:
        vtk = VTKOutput(ma=mesh,
                        coefs=[gfu, flux],
                        names = ["temperature", "flux"],
                        filename="/tmp/thermo_result",
                        subdivision=3)
        vtk.Do()
    
    result = []

    for x_point in x_values:
        point = mesh(x_point, y_value) 
        value_gfu = gfu(point)
        result.append(value_gfu)

    return result


def generate_data(tumor_params: List[TumorParams], debug_index:int = -1, visualize:bool = False, 
                  write_to_csv:bool = True):

    count = 100
    x_values = [x * PARAM_C / count for x in range(count + 1)]
    y_value = PARAM_A

    data = []

    if (debug_index < 0):
        for param in tumor_params:
            result = calculate(param, x_values, y_value, visualize)
            data.append(result)
    else:
        calculate(tumor_params[debug_index], x_values, y_value, visualize)

    if write_to_csv and debug_index < 0:
        with open("/home/user/workspace/project/docs/dissertation/data/thermo_2d_border_data.csv", "w") as data_csv:
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
