# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(settings, model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("expected input shall be a model object")

    return ComputeCFLProcess(model, settings["Parameters"])


class ComputeCFLProcess(Kratos.Process):
    def __init__(self, model, params):
        Kratos.Process.__init__(self)
        default_settings = Kratos.Parameters("""
            {
                "model_part_name"      : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
                "echo_level"           : 0
            }
            """)
        params.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[params["model_part_name"].GetString()]
        self.echo_level = params["echo_level"].GetInt()

    def Execute(self):
        KratosCFD.FluidCharacteristicNumbersUtilities.CalculateLocalCFL(self.model_part)
        if self.echo_level > 0:
            Kratos.Logger.PrintInfo("CFLUtility", "Calculated CFL numebers of elements in {:s}.".format(self.model_part.FullName()))


