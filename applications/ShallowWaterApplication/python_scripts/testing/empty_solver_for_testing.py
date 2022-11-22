# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.ShallowWaterApplication as SW

def CreateSolver(model, custom_settings):
    return EmptySolverForTesting(model, custom_settings)

class EmptySolverForTesting(PythonSolver):
    def __init__(self, model, settings):
        """A solver with the minimal methods to run an analysis."""
        super().__init__(model, settings)
        self.model_part = self.model.CreateModelPart(self.settings["model_part_name"].GetString())
        self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.settings["domain_size"].GetInt())
        self.model_part.ProcessInfo.SetValue(KM.GRAVITY_Z, self.settings["gravity"].GetDouble())
        self.EstimateDeltaTimeUtility = SW.EstimateTimeStepUtility(self.GetComputingModelPart(), self.settings["time_stepping"])

    def ImportModelPart(self):
        self._ImportModelPart(self.model_part,self.settings["model_import_settings"])

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
        self.model_part.AddNodalSolutionStepVariable(SW.FREE_SURFACE_ELEVATION)
        self.model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(SW.TOPOGRAPHY)
        self.model_part.AddNodalSolutionStepVariable(SW.BATHYMETRY)
        self.model_part.AddNodalSolutionStepVariable(SW.MANNING)

    def AddDofs(self):
        formulation_variables = self.settings["formulation_variables"].GetString()
        if formulation_variables == "conservative":
            KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.model_part)
            KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.model_part)
            KM.VariableUtils().AddDof(SW.HEIGHT, self.model_part)
        elif formulation_variables == "primitive":
            KM.VariableUtils().AddDof(KM.VELOCITY_X, self.model_part)
            KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.model_part)
            KM.VariableUtils().AddDof(SW.HEIGHT, self.model_part)
        elif formulation_variables == "boussinesq":
            KM.VariableUtils().AddDof(KM.VELOCITY_X, self.model_part)
            KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.model_part)
            KM.VariableUtils().AddDof(SW.HEIGHT, self.model_part)

    def GetMinimumBufferSize(self):
        return 2

    def GetComputingModelPart(self):
        return self.model_part

    def AdvanceInTime(self, current_time):
        new_time = current_time + self.EstimateDeltaTimeUtility.Execute()
        self.model_part.CloneTimeStep(new_time)
        self.model_part.ProcessInfo[KM.STEP] += 1
        return new_time

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "solver_type"              : "empty_solver_for_testing",
            "model_part_name"          : "model_part",
            "domain_size"              : 2,
            "gravity"                  : 9.81,
            "formulation_variables"    : "conservative",
            "model_import_settings"    : {
                "input_type"               : "use_input_model_part"
            },
            "time_stepping"            : {
                "automatic_time_step"      : false,
                "time_step"                : 0.01
            }
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings
