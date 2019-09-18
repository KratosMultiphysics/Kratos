import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication


def Factory(settings, Model):
    if (type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (type(Model) != KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    turbulence_models_list = ["k_epsilon"]

    model_type = settings["model_type"].GetString()
    if not model_type in turbulence_models_list:
        msg = "Uknown turbulence \"model_type\" name : \"" + model_type
        msg += "\".\nSupported \"model_type\" names: " + turbulence_models_list.__str__(
        )
        raise Exception(msg)

    if model_type == "k_epsilon":
        from KratosMultiphysics.RANSModellingApplication.evm_k_epsilon_configuration import TurbulenceKEpsilonConfiguration
        return TurbulenceKEpsilonConfiguration(Model, settings)