import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeFroudeProcess(Model, settings["Parameters"])

class ComputeFroudeProcess(KM.Process):
    """Compute the Froude number before the output step."""

    def __init__(self, Model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
        {
            "model_part_name"      : "please_specify_model_part_name",
            "save_as_historical"   : false,
            "dry_height_threshold" : 1e-3
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.save_as_historical = settings["save_as_historical"].GetBool()
        self.dry_height_threshold = settings["dry_height_threshold"].GetDouble()

    def ExecuteInitialize(self):
        if not self.save_as_historical:
            KM.VariableUtils().SetNonHistoricalVariableToZero(SW.FROUDE, self.model_part.Nodes)

    def ExecuteBeforeOutputStep(self):
        if self.save_as_historical:
            SW.ShallowWaterUtilities().ComputeFroude(self.model_part, self.dry_height_threshold)
        else:
            SW.ShallowWaterUtilities().ComputeFroudeNonHistorical(self.model_part, self.dry_height_threshold)
