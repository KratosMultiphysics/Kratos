# importing the Kratos Library
from . import structural_response
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication.responses.shape.linear_response as slr
import KratosMultiphysics.OptimizationApplication.responses.shape.plane_symmetry_response as spy
import KratosMultiphysics.OptimizationApplication.responses.partitioning_responses as prr

def CreateResponseFunction(response_name,response_type,response_settings,model):

    if response_type == "mass": 
        return structural_response.MassResponseFunction(response_name,response_settings,model)
    elif response_type == "linear":
        return slr.LinearResponseFunction(response_name,response_settings,model)
    elif response_type == "plane_symmetry":
        return spy.PlaneSymmetryResponseFunction(response_name,response_settings,model)
    elif response_type == "interface":
        return prr.InterfaceResponseFunction(response_name,response_settings,model)      
    elif response_type == "partition_mass":
        return prr.PartitionMassResponseFunction(response_name,response_settings,model)                     