import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication

def Factory(model, settings):
    turbulence_models_list = ["k_epsilon"]

    model_type = settings["model_type"].GetString()
    if not model_type in turbulence_models_list:
        msg = "Uknown turbulence model_type " + model_type
        msg += ". Supported model_types: " + turbulence_models_list
        raise Exception(msg)

    if model_type == "k_epsilon":
        from evm_k_epsilon_configuration import TurbulenceKEpsilonConfiguration
        return TurbulenceKEpsilonConfiguration(model, settings)