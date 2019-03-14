import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as Shallow

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveGeneratorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class WaveGeneratorProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"   : "model_part",
            "formulation"       : "reduced_variables or conserved_variables",
            "wave_length"       : 10.0,
            "wave_period"       : 10.0,
            "wave_height"       : 1.0
        }
        """)

        # Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[settings["model_part_name"].GetString()]

        # Creation of the parameters for the c++ process
        free_surface_parameters = KratosMultiphysics.Parameters("""{}""")
        free_surface_parameters.AddValue("amplitude", settings["wave_height"])
        free_surface_parameters.AddValue("period", settings["wave_period"])
        free_surface_parameters.AddEmptyValue("phase_shift").SetDouble(0.0)
        free_surface_parameters.AddEmptyValue("vertical_shift").SetDouble(0.0)

        self.free_surface_process = Shallow.ApplySinusoidalFunctionToScalar(self.model_part, Shallow.FREE_SURFACE_ELEVATION, free_surface_parameters)
        self.variables_utility = Shallow.ShallowWaterVariablesUtility(self.model_part)

    def ExecuteInitializeSolutionStep(self):

        self.free_surface_process.ExecuteInitializeSolutionStep()
        self.variables_utility.ComputeHeightFromFreeSurface()
        KratosMultiphysics.VariableUtils().ApplyFixity(Shallow.HEIGHT, True, self.model_part.Nodes)
