from typing import Any

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
import KratosMultiphysics.OptimizationApplication.filtering.python_solvers_wrapper_implicit_filters as implicit_filter_solvers
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

# Importing the base class
from KratosMultiphysics.analysis_stage_with_solver import AnalysisStageWithSolver

class HelmholtzAnalysis(AnalysisStageWithSolver):
    """
    This class is the main-script of the implicit filtering analysis
    It can be imported and used as "black-box"
    """
    def __init__(self, model: KM.Model, project_parameters: KM.Parameters):
        super().__init__(model, project_parameters)
        self.__source_data: ContainerExpressionTypes = None
        self.__neighbour_entities: 'dict[Any, KM.Expression.NodalExpression]' = {}

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not already in the model) """
        return implicit_filter_solvers.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[Helmholtz/Sobolev Analysis]:: "

    def _GetComputingModelPart(self):
        return self._GetSolver().GetComputingModelPart()

    def _SetSolverMode(self, invers_mode:bool = False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.COMPUTE_HELMHOLTZ_INVERSE, invers_mode)

    def _SetHelmHoltzSourceMode(self, integrated_field=False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_INTEGRATED_FIELD, integrated_field)

    #### Public user interface functions ####
    def Initialize(self):
        super().Initialize()
        self.InitializeFilterModelPart()

    def InitializeFilterModelPart(self):
        self._SetSolverMode()
        self._GetSolver().SetFilterRadius(self._GetSolver().GetFilterRadius())
        self._SetHelmHoltzSourceMode()

    def RunSolver(self):
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()

    def FilterField(self, unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        self.__AssignDataExpressionToNodalSource(unfiltered_field)
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()
        return self.__AssignNodalSolutionToDataExpression()

    def FilterIntegratedField(self, unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        self.__AssignDataExpressionToNodalSource(unfiltered_field)
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(True)
        self.RunSolver()
        return self.__AssignNodalSolutionToDataExpression()

    def UnFilterField(self, filtered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        self.__AssignDataExpressionToNodalSource(filtered_field)
        self._SetSolverMode(True)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()
        return self.__AssignNodalSolutionToDataExpression()

    def AssignExpressionDataToNodalSolution(self, data_exp: ContainerExpressionTypes) -> None:
        mapped_values = KM.Expression.NodalExpression(data_exp.GetModelPart())
        if isinstance(data_exp, KM.Expression.NodalExpression):
            mapped_values = data_exp
        else:
            KOA.ExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, data_exp)
        KM.Expression.VariableExpressionIO.Write(mapped_values, self._GetSolver().GetSolvingVariable(), True)

    def __AssignDataExpressionToNodalSource(self, data_exp: ContainerExpressionTypes):
        self.__source_data = data_exp

        # it is better to work on the model part of the data_exp rather than the internal
        # model part created with ConnectivityPreserveModelPart because, then all the outputs will
        # be written to one vtu file automatically, otherwise, there will be outputs distributed among
        # the helmholtz model part and the original model part. This is safer because the helmholtz modelpart
        # is created using the ConnectivityPreserveModeller which preserves the same nodes and condition/element
        # data containers.
        mapped_values = KM.Expression.NodalExpression(data_exp.GetModelPart())
        if isinstance(data_exp, KM.Expression.NodalExpression):
            mapped_values = data_exp
        else:
            KOA.ExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, data_exp)

        KM.Expression.VariableExpressionIO.Write(mapped_values, self._GetSolver().GetSourceVariable(), False)

    def __AssignNodalSolutionToDataExpression(self) -> ContainerExpressionTypes:
        if self.__source_data is None:
            raise RuntimeError("The __AssignDataExpressionToNodalSource should be called first.")

        # it is better to work on the model part of the data_exp rather than the internal
        # model part created with ConnectivityPreserveModelPart because, then all the outputs will
        # be written to one vtu file automatically, otherwise, there will be outputs distributed among
        # the helmholtz model part and the original model part. This is safer because the helmholtz modelpart
        # is created using the ConnectivityPreserveModeller which preserves the same nodes and condition/element
        # data containers.
        nodal_solution_field = KM.Expression.NodalExpression(self.__source_data.GetModelPart())

        KM.Expression.VariableExpressionIO.Read(nodal_solution_field, self._GetSolver().GetSolvingVariable(), True)

        if isinstance(self.__source_data, KM.Expression.NodalExpression):
            return nodal_solution_field.Clone()
        else:
            mapped_entity_solution_field = self.__source_data.Clone()
            KOA.ExpressionUtils.MapNodalVariableToContainerVariable(mapped_entity_solution_field, nodal_solution_field, self.__GetNeighbourEntities(self.__source_data))
            return mapped_entity_solution_field

    def __GetNeighbourEntities(self, data_exp: ContainerExpressionTypes) -> ContainerExpressionTypes:
        # following makes the number of neighbours computation to be executed once
        # per given container, hence if the mesh element/connectivity changes
        # this computation needs to be redone. Especially in the case if MMG is
        # used for re-meshing.
        key = data_exp.GetContainer()
        if key not in  self.__neighbour_entities.keys():
            self.__neighbour_entities[key] = KM.Expression.NodalExpression(data_exp.GetModelPart())
            if isinstance(data_exp, KM.Expression.ElementExpression):
                KOA.ExpressionUtils.ComputeNumberOfNeighbourElements(self.__neighbour_entities[key])
            else:
                KOA.ExpressionUtils.ComputeNumberOfNeighbourConditions(self.__neighbour_entities[key])
        return self.__neighbour_entities[key]