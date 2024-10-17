# importing the Kratos Library
from . import structural_responses
from . import additive_manufacturing_responses

def CreateResponseFunction(response_name,response_type,response_settings, response_analysis,model):

    if response_type == "strain_energy": 
        return structural_responses.StrainEnergyResponseFunction(response_name,response_settings,response_analysis,model)

    elif response_type == "stress":
        return structural_responses.StressResponseFunction(response_name,response_settings,response_analysis,model)

    elif response_type == "partition_interface_stress":
        return additive_manufacturing_responses.PartitionInterfaceStressResponseFunction(response_name,response_settings,model)        
