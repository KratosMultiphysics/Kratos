# importing the Kratos Library
from . import structural_response
from . import partitioning_responses

def CreateResponseFunction(response_name,response_type,response_settings, response_analysis,model):

    if response_type == "strain_energy": 
        return structural_response.StrainEnergyResponseFunction(response_name,response_settings,response_analysis,model)

    elif response_type == "stress":
        return structural_response.StressResponseFunction(response_name,response_settings,response_analysis,model)

    elif response_type == "partition_interface_stress":
        return partitioning_responses.PartitionInterfaceStressResponseFunction(response_name,response_settings,model)        

    # elif response_type == "mass":
    #     return structural_response.MassResponseFunction(response_id, response_settings, model)

    # elif response_type == "eigenfrequency":
    #     return structural_response.EigenFrequencyResponseFunction(response_id, response_settings, model)

    # elif response_type == "adjoint_nodal_displacement":
    #     return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    # elif response_type == "adjoint_linear_strain_energy":
    #     return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    # elif response_type == "adjoint_local_stress":
    #     return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    # elif response_type == "adjoint_max_stress":
    #     return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    # elif response_type == "adjoint_nodal_reaction":
    #     return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    # else:
    #     raise NameError("The type of the following response function is not specified: "+ response_id +
    #                     ".\nAvailable types are: 'mass', 'strain_energy', 'eigenfrequency', 'adjoint_nodal_displacement', 'adjoint_linear_strain_energy', 'adjoint_local_stress', 'adjoint_max_stress', 'adjoint_nodal_reaction'." )
