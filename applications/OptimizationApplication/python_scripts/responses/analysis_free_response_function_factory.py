# importing the Kratos Library
from . import structural_responses
from . import additive_manufacturing_responses

def CreateResponseFunction(response_name,response_type,response_settings,model):

    if response_type == "mass":
        return structural_responses.MassResponseFunction(response_name,response_settings,model)
    elif response_type == "interface":
        return additive_manufacturing_responses.InterfaceResponseFunction(response_name,response_settings,model)
    elif response_type == "partition_mass":
        return additive_manufacturing_responses.PartitionMassResponseFunction(response_name,response_settings,model)
    elif response_type == "am_max_overhang_angle":
        return additive_manufacturing_responses.MaxOverhangAngleResponseFunction(response_name,response_settings,model)