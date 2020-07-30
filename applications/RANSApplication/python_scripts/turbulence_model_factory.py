import KratosMultiphysics
import KratosMultiphysics.RANSApplication

from KratosMultiphysics.RANSApplication.evm_k_epsilon_configuration import TurbulenceKEpsilonConfiguration

def Factory(settings, Model):
    if (not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (not isinstance(Model, KratosMultiphysics.Model)):
        raise Exception("expected input shall be a Model object")
    turbulence_models_list = ["k_epsilon"]

    model_type = settings["model_type"].GetString()
    if not model_type in turbulence_models_list:
        msg = "Unknown turbulence \"model_type\" name : \"" + model_type
        msg += "\".\nSupported \"model_type\" names:"
        msg += "\n\t".join(turbulence_models_list)
        raise Exception(msg)

    if model_type == "k_epsilon":
        turbulence_model = TurbulenceKEpsilonConfiguration(Model, settings)
        return turbulence_model
