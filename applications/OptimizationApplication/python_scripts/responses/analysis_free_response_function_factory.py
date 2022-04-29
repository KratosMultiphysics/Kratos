# importing the Kratos Library
from . import structural_response
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication.response_functions.response_function_factory as sho_response_factory
import KratosMultiphysics.OptimizationApplication.responses.shape.linear_response as slr
import KratosMultiphysics.OptimizationApplication.responses.shape.plane_symmetry_response as spy


def CreateResponseFunction(response_name,response_type,response_settings,model):

    if response_type == "mass": 
        return structural_response.MassResponseFunction(response_name,response_settings,model)
    elif response_type == "linear":
        return slr.LinearResponseFunction(response_name,response_settings,model)
    elif response_type == "plane_symmetry":
        return spy.PlaneSymmetryResponseFunction(response_name,response_settings,model)        