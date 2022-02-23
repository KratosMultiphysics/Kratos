# importing the Kratos Library
from . import structural_response

def CreateResponseFunction(response_type,response_settings, response_analyzer,response_analyzer_model_part,model):

    structural_primal_analysis_required_responses = ["strain_energy","eigenfrequency","adjoint_nodal_displacement","adjoint_linear_strain_energy",
                                          "adjoint_local_stress","adjoint_max_stress","adjoint_nodal_reaction"]

    if response_type in structural_primal_analysis_required_responses and (response_analyzer is None or type(response_analyzer) is not structural_response.StructuralMechanicsAnalysis):
        raise RuntimeError("CreateResponseFunction: Response {} requires analysis type of StructuralMechanicsAnalysis which is not defined !".format(response_type))  

    if response_type == "strain_energy": 
        return structural_response.StrainEnergyResponseFunction(response_settings,response_analyzer,response_analyzer_model_part,model)

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
