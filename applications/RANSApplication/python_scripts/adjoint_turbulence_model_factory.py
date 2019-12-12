import KratosMultiphysics
import KratosMultiphysics.RANSApplication

from KratosMultiphysics.RANSApplication.adjoint_evm_k_epsilon_configuration import AdjointTurbulenceKEpsilonConfiguration


def Factory(settings, Model):
    if (not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (not isinstance(Model, KratosMultiphysics.Model)):
        raise Exception("expected input shall be a Model object")
    turbulence_models_list = ["k_epsilon"]

    model_type = settings["model_type"].GetString()
    if model_type not in turbulence_models_list:
        msg = "Uknown adjoint turbulence \"model_type\" name : \"" + model_type
        msg += "\".\nSupported \"model_type\" names: " + turbulence_models_list.__str__(
        )
        raise Exception(msg)

    if model_type == "k_epsilon":
        return AdjointTurbulenceKEpsilonConfiguration(Model, settings)