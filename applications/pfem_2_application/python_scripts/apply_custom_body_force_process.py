import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomBodyForceProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"     : "please_specify_model_part_name",
                "variable_name"       : "BODY_FORCE",
                "body_force_type"     : "vortex",
                "Parameters"          : {}
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())

        if settings["body_force_type"].GetString() == "vortex":
            from manufactured_vortex import ManufacturedVortex
            self.body_force = ManufacturedVortex(settings["Parameters"])
        else:
            msg = "Body force type not implemented"
            raise Exception(msg)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        self.CalculateBodyForce()

    def CalculateBodyForce(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        for node in self.model_part.Nodes:
            value = self.body_force.BodyForce(node.X, node.Y, node.Z, current_time)
            node.SetSolutionStepValue(self.variable, value)
