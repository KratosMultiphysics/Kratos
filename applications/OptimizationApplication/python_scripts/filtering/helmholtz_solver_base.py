import abc
import typing

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

# Other imports
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

class HelmholtzSolverBase(PythonSolver):
    """The base class for Helmholtz-based solvers.

    This class defines the user interface to Helmholtz solvers.

    """
    def __init__(self, model: KM.Model, custom_settings: KM.Parameters):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super().__init__(model, custom_settings)

        # Either retrieve the model part from the model or create a new one
        self.__filtering_model_part_name = self.settings["model_part_name"].GetString()

        if self.__filtering_model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        # the __filtering_model_part_name can be a sub model part. Hence we need the main model part
        # for variable and dof addition.
        root_model_part_name = self.__filtering_model_part_name.split(".")[0]
        if self.model.HasModelPart(root_model_part_name):
            self.origin_root_model_part = self.model.GetModelPart(root_model_part_name)
        else:
            self.origin_root_model_part = model.CreateModelPart(root_model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')
        self.origin_root_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)

        # we don't create the helmholtz model part here, because
        # we need to have unique model parts for different filtering variable types
        # hence, it cannot be determined at this point.
        self.__helmholtz_model_part: 'typing.Optional[KM.ModelPart]' = None

        # Get the filter radius
        self.filter_radius = self.settings["filter_radius"].GetDouble()
        self.filter_type = self.settings["filter_type"].GetString()

        KM.Logger.PrintInfo("::[HelmholtzSolverBase]:: Construction finished")

    @classmethod
    def GetDefaultParameters(cls) -> KM.Parameters:
        this_defaults = KM.Parameters("""{
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

    def AddVariables(self) -> None:
        self.GetOriginRootModelPart().AddNodalSolutionStepVariable(self._GetSolvingVariable())
        KM.Logger.PrintInfo("::[HelmholtzSolverBase]:: Variables ADDED.")

    def AddDofs(self) -> None:
        if isinstance(self._GetSolvingVariable(), KM.DoubleVariable):
            KM.VariableUtils().AddDof(self._GetSolvingVariable(), self.GetOriginRootModelPart())
        elif isinstance(self._GetSolvingVariable(), KM.Array1DVariable3):
            KM.VariableUtils().AddDof(KM.KratosGlobals.GetVariable(f"{self._GetSolvingVariable().Name()}_X"), self.GetOriginRootModelPart())
            KM.VariableUtils().AddDof(KM.KratosGlobals.GetVariable(f"{self._GetSolvingVariable().Name()}_Y"), self.GetOriginRootModelPart())
            KM.VariableUtils().AddDof(KM.KratosGlobals.GetVariable(f"{self._GetSolvingVariable().Name()}_Z"), self.GetOriginRootModelPart())
        else:
            raise RuntimeError("Unsupported solving variable type.")
        KM.Logger.PrintInfo("::[HelmholtzSolverBase]:: DOFs ADDED.")

    def AdvanceInTime(self, current_time: float) -> float:
        dt = self.settings["time_stepping"]["time_step"].GetDouble()
        new_time = current_time + dt
        self.__helmholtz_model_part.ProcessInfo[KM.STEP] += 1
        self.__helmholtz_model_part.CloneTimeStep(new_time)

        return new_time

    def Initialize(self) -> None:
        self._GetSolutionStrategy().Initialize()
        neighbours_exp = KM.Expression.NodalExpression(self.__helmholtz_model_part)
        KOA.ExpressionUtils.ComputeNumberOfNeighbourElements(neighbours_exp)
        KM.Expression.VariableExpressionIO.Write(neighbours_exp, KM.NUMBER_OF_NEIGHBOUR_ELEMENTS, False)
        KM.Logger.PrintInfo("::[HelmholtzSolverBase]:: Finished initialization.")

    def InitializeSolutionStep(self) -> None:
        self._GetSolutionStrategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self) -> None:
        self._GetSolutionStrategy().FinalizeSolutionStep()

    def Predict(self) -> None:
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self) -> None:
        is_converged = bool(self._GetSolutionStrategy().Solve())
        return is_converged

    def SetEchoLevel(self, level: int) -> None:
        self._GetSolutionStrategy().SetEchoLevel(level)

    def GetEchoLevel(self) -> int:
        self._GetSolutionStrategy().GetEchoLevel()

    def Clear(self) -> None:
        self._GetSolutionStrategy().Clear()

    def ImportModelPart(self) -> None:
        self._ImportModelPart(self.origin_root_model_part, self.settings["model_import_settings"])

    def GetComputingModelPart(self) -> KM.ModelPart:
        return self.__helmholtz_model_part

    def GetOriginRootModelPart(self) -> KM.ModelPart:
        return self.origin_root_model_part

    def GetOriginModelPart(self) -> KM.ModelPart:
        return self.model[self.__filtering_model_part_name]

    def GetFilterType(self) -> str:
        return self.filter_type

    def GetFilterRadius(self) -> str:
        return self.filter_radius

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
        return KM.ResidualBasedBlockBuilderAndSolver(linear_solver)

    def _CreateLinearSolver(self):
        return ConstructSolver(self.settings["linear_solver_settings"])

    def _CreateScheme(self):
        return KM.ResidualBasedIncrementalUpdateStaticScheme()

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        return KM.ResidualBasedLinearStrategy(computing_model_part,
                                                              scheme,
                                                              builder_and_solver,
                                                              False,
                                                              False,
                                                              False,
                                                              False)

    def _GetContainerTypeNumNodes(self, container: 'typing.Union[KM.ConditionsArray, KM.ElementsArray]') -> int:
        num_nodes = None
        for cont_type in container:
            return len(cont_type.GetNodes())
        return num_nodes

    def _IsSurfaceContainer(self, container: 'typing.Union[KM.ConditionsArray, KM.ElementsArray]') -> bool:
        return any(map(lambda x: x.GetGeometry().WorkingSpaceDimension() != x.GetGeometry().LocalSpaceDimension(), container))

    def GetComputingModelPartName(self) -> str:
        model_part_name = self.__filtering_model_part_name.replace(".", "_")
        for container in self.GetContainers():
            if isinstance(container, KM.ConditionsArray):
                model_part_name += "_conditions"
            elif isinstance(container, KM.ElementsArray):
                model_part_name += "_elements"
            else:
                raise RuntimeError("Unsupported container type provided.")
        return model_part_name + f"_{self.GetConditionName()}_{self.GetElementName()}"

    @abc.abstractmethod
    def _GetComputingModelPartName(self) -> str:
        pass

    @abc.abstractmethod
    def _FillComputingModelPart(self) -> None:
        pass

    @abc.abstractmethod
    def _GetSolvingVariable(self) -> SupportedSensitivityFieldVariableTypes:
        pass

    def _GetSourceVariable(self) -> SupportedSensitivityFieldVariableTypes:
        if isinstance(self._GetSolvingVariable(), KM.DoubleVariable):
            return KOA.HELMHOLTZ_SCALAR_SOURCE
        elif isinstance(self._GetSolvingVariable(), KM.Array1DVariable3):
            return KOA.HELMHOLTZ_VECTOR_SOURCE
        else:
            raise RuntimeError("Unsupported solving variable.")

    def PrepareModelPart(self) -> str:
        computing_model_part_name = self._GetComputingModelPartName()

        if self.model.HasModelPart(computing_model_part_name):
            # a same model part with same entities were already created by some other filter
            # hence reusing it
            self.__helmholtz_model_part = self.model[computing_model_part_name]

            # we increase the number of filters using the model part
            # this is required to know whether we need to apply boundary conditions
            # every time filter is used so that we don't mix-up different BCs from
            # different filters
            self.__helmholtz_model_part[KOA.NUMBER_OF_HELMHOLTZ_FILTERS] += 1
        else:
            # the model part needs to be created.
            self.__helmholtz_model_part = self.model.CreateModelPart(computing_model_part_name)
            KOA.OptimizationUtils.CopySolutionStepVariablesList(self.__helmholtz_model_part, self.GetOriginRootModelPart())
            self.__helmholtz_model_part[KOA.NUMBER_OF_HELMHOLTZ_FILTERS] = 1

            # now fill with the appropriate elements and conditions
            self._FillComputingModelPart()

