# importing the Kratos Library
from turtle import distance
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.wave_solver import WaveSolver

def CreateSolver(model, custom_settings):
    return EmbeddedPrimitiveSolver(model, custom_settings)


class EmbeddedPrimitiveSolver(WaveSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.element_name, self.condition_name, self.min_buffer_size = self._GetFormulationSettings()

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA) # We need this for the parallel redistance

    def Initialize(self):
        super().Initialize()

        # Initialize the distance field. NOTE: the initial condition must be set at Process.ExecuteInitialize
        KM.VariableUtils().CopyVariable(SW.HEIGHT, KM.DISTANCE, self.main_model_part.Nodes)

        # Instantiate the distance convection process
        self.GetDistanceConvectionProcess()

        # Call the redistancing process to make sure the initial condition is a distance field
        self.GetDistanceReinitializationProcess().Execute()

        # Set the distance modification process
        self.GetDistanceModificationProcess().ExecuteInitialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        # velocity = self.GetComputingModelPart().GetNode(88).GetSolutionStepValue(KM.VELOCITY)
        # KM.VariableUtils().SetVariable(KM.VELOCITY_X, velocity[0], self.GetComputingModelPart().Nodes)
        # KM.VariableUtils().SetVariable(KM.VELOCITY_Y, 0.0, self.GetComputingModelPart().Nodes)
        for node in self.GetComputingModelPart().Nodes:
            node.Fix(KM.VELOCITY_Y)

        # Perform distance convection
        self.GetDistanceConvectionProcess().Execute()

        # Recompute the distance field according to the new level-set position
        self.GetDistanceReinitializationProcess().Execute()

        # Do the distance correction and set the fixity in the "internal" nodes
        self.GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # TODO: Here we must implement the FM-ALE (or similar) type operations

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # Restore the fixity to its original status
        self.GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

    def _GetFormulationSettings(self):
        scheme = self.settings["time_integration_scheme"].GetString()
        order = self.settings["time_integration_order"].GetInt()
        if scheme == "bdf":
            element_name = "EmbeddedPrimitiveElement"
            condition_name = "WaveCondition"
            buffer_size = order + 1
        else:
            raise Exception('The possible "time_integration_scheme" are "bdf" and "crank_nicolson"')
        return element_name, condition_name, buffer_size

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""{
            "distance_modification_settings" : {},
            "distance_convection_settings" : {},
            "distance_reinitialization_type" : "variational"
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process

    def GetDistanceReinitializationProcess(self):
        if not hasattr(self, '_distance_reinitialization_process'):
            self._distance_reinitialization_process = self._CreateDistanceReinitializationProcess()
        return self._distance_reinitialization_process

    def GetDistanceConvectionProcess(self):
        if not hasattr(self, '_level_set_convection_process'):
            self._level_set_convection_process = self._CreateDistanceConvectionProcess()
        return self._level_set_convection_process

    @classmethod
    def __GetDistanceModificationDefaultSettings(self):
        return KM.Parameters(r'''{
            "model_part_name": "",
            "distance_threshold": 1e-3,
            "continuous_distance": true,
            "check_at_each_time_step": true,
            "avoid_almost_empty_elements": true,
            "deactivate_full_negative_elements": true,
            "full_negative_elements_fixed_variables_list" : ["HEIGHT","VELOCITY"]
        }''')

    def __CreateDistanceModificationProcess(self):
        # Set the distance modification settings according to the level set type
        # Note that the distance modification process is applied to the volume model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(
            self.__GetDistanceModificationDefaultSettings())
        volume_part_name = self.settings["model_part_name"].GetString()
        distance_modification_settings["model_part_name"].SetString(volume_part_name)
        return KratosFluid.DistanceModificationProcess(self.model, distance_modification_settings)

    def _CreateDistanceConvectionProcess(self):
        # Construct the level set convection process
        computing_model_part = self.GetComputingModelPart()
        linear_solver = self._GetLevelsetLinearSolver()
        distance_convection_settings = self.settings["distance_convection_settings"]
        level_set_convection_process = KM.LevelSetConvectionProcess2D(
            computing_model_part,
            linear_solver,
            distance_convection_settings)
        return level_set_convection_process

    def _GetLevelsetLinearSolver(self):
        # A linear solver configured specifically for the level-set convection process
        if not hasattr(self, '_levelset_linear_solver'):
            # TODO: add customized configuration
            self._levelset_linear_solver = self._CreateLinearSolver()
        return self._levelset_linear_solver

    def _CreateDistanceReinitializationProcess(self):
        reinitialization_type = self.settings["distance_reinitialization_type"].GetString()
        if reinitialization_type == "variational":
            maximum_iterations = 2  # TODO: Make this user-definable
            linear_solver = self._GetLevelsetLinearSolver()
            computing_model_part = self.GetComputingModelPart()
            distance_reinitialization_process = KM.VariationalDistanceCalculationProcess2D(
                computing_model_part,
                linear_solver,
                maximum_iterations,
                KM.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        elif reinitialization_type == "parallel":
            parallel_distance_settings = KM.Parameters("""{
                "max_levels" : 25,
                "max_distance" : 1.0,
                "calculate_exact_distances_to_plane" : true
            }""")
            distance_reinitialization_process = KM.ParallelDistanceCalculationProcess2D(
                self.main_model_part,
                parallel_distance_settings)

        elif reinitialization_type == "none":
            distance_reinitialization_process = KM.Process()
        else:
            raise Exception("Please use a valid distance reinitialization type or set it as \'none\'. Valid types are: \'variational\' and \'parallel\'.")

        return distance_reinitialization_process
