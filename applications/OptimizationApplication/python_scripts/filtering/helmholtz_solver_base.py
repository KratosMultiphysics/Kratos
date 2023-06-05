# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.OptimizationApplication as KOA

# Other imports
from KratosMultiphysics.python_solver import PythonSolver

class HelmholtzSolverBase(PythonSolver):
    """The base class for Helmholtz-based solvers.

    This class defines the user interface to Helmholtz solvers.

    """
    def __init__(self, model: KratosMultiphysics.Model, custom_settings: KratosMultiphysics.Parameters):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super().__init__(model, custom_settings)

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.original_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.original_model_part = model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.original_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Create Helmholtz model part
        self.helmholtz_model_part = self.model.CreateModelPart(self.original_model_part.Name+"_helmholtz_filter_mdp")

        # Get the filter radius
        self.filter_radius = self.settings["filter_radius"].GetDouble()

        KratosMultiphysics.Logger.PrintInfo("::[HelmholtzSolverBase]:: Construction finished")

    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type"           : "helmholtz_solver_base",
            "domain_size"           : -1,
            "filter_type"     : "",
            "filter_radius"     : 0.0,
            "model_part_name"       : "",
            "time_stepping" : {
                "time_step"       : 1.0
            },
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "linear_solver_settings" : {
                "solver_type" : "amgcl",
                "smoother_type":"ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 5000
            },
            "material_properties": {}
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    #### Public user interface functions ####

    def AdvanceInTime(self, current_time) -> float:
        dt = self.settings["time_stepping"]["time_step"].GetDouble()
        new_time = current_time + dt
        self.helmholtz_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.helmholtz_model_part.CloneTimeStep(new_time)

        return new_time

    def Initialize(self):
        self._GetSolutionStrategy().Initialize()
        KOA.ImplicitFilterUtils.CalculateNodeNeighbourCount(self.helmholtz_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[HelmholtzSolverBase]:: Finished initialization.")

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        is_converged = bool(self._GetSolutionStrategy().Solve())
        return is_converged

    def SetEchoLevel(self, level):
        self._GetSolutionStrategy().SetEchoLevel(level)

    def GetEchoLevel(self):
        self._GetSolutionStrategy().GetEchoLevel()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.original_model_part, self.settings["model_import_settings"])

    def GetComputingModelPart(self):
        return self.helmholtz_model_part

    #### Protected functions ####

    #### Specific internal functions ####

    def _GetScheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._CreateScheme()
        return self._solution_scheme

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        return KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)

    def _CreateLinearSolver(self):
        from KratosMultiphysics.python_linear_solver_factory import ConstructSolver
        return ConstructSolver(self.settings["linear_solver_settings"])

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              scheme,
                                                              builder_and_solver,
                                                              False,
                                                              False,
                                                              False,
                                                              False)

    def _AssignProperties(self, parameters: KratosMultiphysics.Parameters):
        KOA.ImplicitFilterUtils.AssignProperties(self.GetComputingModelPart(), parameters)

