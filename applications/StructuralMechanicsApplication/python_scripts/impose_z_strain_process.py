import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]
    default_params = KratosMultiphysics.Parameters("""{
        "model_part_name" : "please_specify_model_part_name",
        "z_strain_value" : 0.01
    }""")
    process_settings.AddMissingParameters(default_params)

    model_part = Model[process_settings["model_part_name"].GetString()]

    return SMA.ImposeZStrainProcess(model_part, process_settings)