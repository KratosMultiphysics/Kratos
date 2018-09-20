from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import structural_response
import structural_response_global_finite_differencing

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

    elif response_type == "global_finite_differencing":
        return structural_response_global_finite_differencing.GlobalFiniteDifferencingResponseFunction(response_id, response_settings, model)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'mass', 'strain_energy', 'eigenfrequency', 'adjoint_nodal_displacement', 'adjoint_linear_strain_energy', 'adjoint_local_stress'." )