import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyAbsorbingBoundaryProcess(model, settings["Parameters"])

class ApplyAbsorbingBoundaryProcess(KM.Process):
    """ApplyAbsorbingBoundaryProcess

    This process sets the DISTANCE variable from every
    node to the nearest boundary condition
    """

    def __init__(self, model, settings):
        """The constructor of the ApplyAbsorbingBoundaryProcess

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The model to be used
        settings -- The ProjectParameters used
        """

        KM.Process.__init__(self)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.processes = []
        self.model_part = model[settings["computing_model_part_name"].GetString()]
        boundaries_names = settings["absorbing_boundaries_list"].GetStringArray()
        for name in boundaries_names:
            boundary_part = self.model_part.GetSubModelPart(name)
            self.processes.append(SW.CalculateDistanceToBoundaryProcess(self.model_part, boundary_part, settings))

        for process in self.processes:
            process.Check()

    def ExecuteInitialize(self):
        KM.VariableUtils().SetVariable(KM.DISTANCE, 1e+38, self.model_part.Nodes)

    def ExecuteBeforeSolutionLoop(self):
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
                "computing_model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "absorbing_boundaries_list" : [],
                "r_squared_treshold"        : 0.99
            }""")
