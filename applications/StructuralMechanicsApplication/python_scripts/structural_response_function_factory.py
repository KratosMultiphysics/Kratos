# importing the Kratos Library
from KratosMultiphysics.StructuralMechanicsApplication import structural_response

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "strain_energy":
        return structural_response.StrainEnergyResponseFunction(response_id, response_settings, model)

    elif response_type == "mass":
        return structural_response.MassResponseFunction(response_id, response_settings, model)

    elif response_type == "eigenfrequency":
        return structural_response.EigenFrequencyResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_nodal_displacement":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_linear_strain_energy":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_local_stress":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_max_stress":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_nodal_reaction":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'mass', 'strain_energy', 'eigenfrequency', 'adjoint_nodal_displacement', 'adjoint_linear_strain_energy', 'adjoint_local_stress', 'adjoint_max_stress', 'adjoint_nodal_reaction'." )
