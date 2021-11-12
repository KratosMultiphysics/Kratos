import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetTopographyProcess(Model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class SetTopographyProcess(KM.Process):

    def __init__(self, Model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "interval"             : [0.0, 1e30],
                "constrained"          : false,
                "value"                : "z",
                "set_mesh_z_to_zero"   : true
            }
            """
            )
        if settings.Has("value"):
            if settings["value"].IsDouble():
                default_settings["value"].SetDouble(0.0)
        settings.ValidateAndAssignDefaults(default_settings)
        settings.AddEmptyValue("variable_name").SetString("TOPOGRAPHY")

        self.model_part = Model[settings["model_part_name"].GetString()]
        if settings["value"].IsString():
            self.depends_on_time = settings["value"].GetString().find('t') != -1
        else:
            self.depends_on_time = False
        self.set_mesh_z_to_zero = settings["set_mesh_z_to_zero"].GetBool()
        process_settings = settings.Clone()
        process_settings.RemoveValue("set_mesh_z_to_zero")

        from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
        self.process = AssignScalarVariableProcess(Model, process_settings)

    def ExecuteInitialize(self):
        self.process.ExecuteInitializeSolutionStep()
        SW.ShallowWaterUtilities().FlipScalarVariable(SW.TOPOGRAPHY, SW.BATHYMETRY, self.model_part)
        if self.set_mesh_z_to_zero:
            SW.ShallowWaterUtilities().SetMeshZCoordinateToZero(self.model_part)
            SW.ShallowWaterUtilities().SetMeshZ0CoordinateToZero(self.model_part)

    def ExecuteInitializeSolutionStep(self):
        if self.depends_on_time:
            self.ExecuteInitialize()
