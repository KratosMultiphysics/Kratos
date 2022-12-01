
from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response

try:
    import KratosMultiphysics.MultilevelMonteCarloApplication
    import KratosMultiphysics.MappingApplication
    import xmc
    import exaqute
    import numpy as np
    from KratosMultiphysics.CompressiblePotentialFlowApplication import stochastic_potential_flow_response
    is_xmc_available = True
except:
    is_xmc_available = False

def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "adjoint_lift_potential_jump":
        return potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
    elif response_type == "stochastic_adjoint_lift_potential_jump":
        if is_xmc_available:
            return stochastic_potential_flow_response.AdjointResponseFunction(response_id, response_settings, model)
        else:
            err_msg = "XMC and its dependencies could not be imported. Please check applications/MultilevelMonteCarloApplication/README.md"
            err_msg += " for installation details. Please make sure that the CompressiblePotentialFlowApplication, MultilevelMonteCarloApplication"
            err_msg += " and the MappingApplication are compiled."
            raise ImportError(err_msg)
    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'adjoint_lift_potential_jump'." )
