# importing the Kratos Library
from . import structural_response
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication.response_functions.response_function_factory as sho_response_factory
from . import linear_response


def CreateResponseFunction(response_name,response_type,response_settings,model):

    if response_type == "mass": 
        return structural_response.MassResponseFunction(response_name,response_settings,model)
    elif response_type == "linear":
        return linear_response.LinearResponseFunction(response_name,response_settings,model)