import sympy
import KratosMultiphysics
from KratosMultiphysics.sympy_fe_utilities import DefineMatrix

from params_dict import params
from shape_functions import DefineShapeFunctionsMatrix
import generate_convective_flux
import generate_diffusive_flux
import generate_source_term
import generate_stabilization_matrix

do_simplifications = False          # Simplify resulting differenctiations
subscales_vector = ["ASGS", "OSS"]  # Subscales types to be computed
dim_vector = [2]                    # Spatial dimensions to be computed
is_explicit = True                  # Explicit or implicit time integration

# template_filename = "" # TODO
# with open(template_filename) as f:
#     outstring = f.read()

for dim in dim_vector:
    # Change dimension accordingly
    params["dim"] = dim

    # Shape functions and Gauss pts. settings
    (n_nodes, n_gauss) = {
        1: (1, 2),
        2: (3, 3),
        3: (8, 8)
    }[dim]

    DN = DefineMatrix('DN', n_nodes, dim)
    mat_N = DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss)