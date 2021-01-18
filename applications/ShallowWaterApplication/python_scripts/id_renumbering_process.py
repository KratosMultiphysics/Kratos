import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return IdRenumberingProcess(model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class IdRenumberingProcess(KM.Process):

    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "renumber_all_model_parts" : true,
                "model_parts_list"         : [],
                "renumber_nodes"           : true,
                "renumber_elements"        : true,
                "renumber_conditions"      : true
            }
            """
            )
        settings.ValidateAndAssignDefaults(default_settings)

        if settings["renumber_all_model_parts"].GetBool():
            self.process = SW.IdRenumberingProcess(model)
        else:
            self.process = SW.IdRenumberingProcess(model, settings["model_parts_list"].GetStringArray())

        self.renumber_nodes = settings["renumber_nodes"].GetBool()
        self.renumber_elements = settings["renumber_elements"].GetBool()
        self.renumber_conditions = settings["renumber_conditions"].GetBool()

    def ExecuteBeforeOutputStep(self):
        if self.renumber_nodes:
            self.process.RenumberNodes()
        if self.renumber_elements:
            self.process.RenumberElements()
        if self.renumber_conditions:
            self.process.RenumberConditions()

    def ExecuteAfterOutputStep(self):
        if self.renumber_nodes:
            self.process.RestoreNodes()
        if self.renumber_elements:
            self.process.RestoreElements()
        if self.renumber_conditions:
            self.process.RestoreConditions()
