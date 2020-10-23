# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.ShallowWaterApplication as SW

def CreateSolver(model, custom_settings):
    return EmptySolverForTesting(model, custom_settings)

class EmptySolverForTesting(PythonSolver):
    def __init__(self, model, settings):
        """ A solver with the minimal methods to run an analysis.

        This class can be used to test the processes.
        """
        self._validate_settings_in_baseclass = True
        super().__init__(model, settings)

        model_part_name = self.settings["model_part_name"].GetString()
        if self.model.HasModelPart(model_part_name):
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            self.model_part = self.model.CreateModelPart(model_part_name)

        self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, self.settings["domain_size"].GetInt())
        self.model_part.ProcessInfo.SetValue(KM.GRAVITY_Z, self.settings["gravity"].GetDouble())

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
        self.model_part.AddNodalSolutionStepVariable(KM.NODAL_H)

    def GetMinimumBufferSize(self):
        return 2

    def GetComputingModelPart(self):
        return self.model_part

    def Initialize(self):
        # The time step utility needs the NODAL_H
        KM.FindNodalHProcess(self.GetComputingModelPart()).Execute()
        self.EstimateDeltaTimeUtility = SW.EstimateDtShallow(self.GetComputingModelPart(), self.settings["time_stepping"])

    def AdvanceInTime(self, current_time):
        dt = self.EstimateDeltaTimeUtility.EstimateDt()
        new_time = current_time + dt
        self.model_part.CloneTimeStep(new_time)
        self.model_part.ProcessInfo[KM.STEP] += 1
        self.model_part.ProcessInfo[KM.DELTA_TIME] = dt
        return new_time

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "solver_type"              : "empty_solver_for_testing",
            "model_part_name"          : "model_part",
            "domain_size"              : 2,
            "gravity"                  : 9.81,
            "model_import_settings"    : {
                "input_type"               : "mdpa",
                "input_filename"           : "unknown_name"
            },
            "time_stepping"            : {
                "automatic_time_step"      : false,
                "time_step"                : 0.01
            }
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings
