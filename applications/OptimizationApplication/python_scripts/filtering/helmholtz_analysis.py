# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
import KratosMultiphysics.OptimizationApplication.filtering.python_solvers_wrapper_implicit_filters as implicit_filter_solvers
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class HelmholtzAnalysis(AnalysisStage):
    """
    This class is the main-script of the implicit filtering analysis
    It can be imported and used as "black-box"
    """
    def __init__(self, model: KM.Model, project_parameters: KM.Parameters):

        self.initialized = False
        super().__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        return implicit_filter_solvers.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[Helmholtz/Sobolev Analysis]:: "

    def _GetComputingModelPart(self):
        return self._GetSolver().GetComputingModelPart()

    def _SetSolverMode(self, invers_mode:bool = False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.COMPUTE_HELMHOLTZ_INVERSE, invers_mode)

    def _SetHelmHoltzSourceMode(self, integrated_field=False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_INTEGRATED_FIELD, integrated_field)

    def _AssignDataExpressionToNodalSource(self, data_exp: ContainerExpressionTypes):
        mapped_values = KM.Expression.NodalNonHistoricalExpression(data_exp.GetModelPart())
        if isinstance(data_exp, KM.Expression.NodalNonHistoricalExpression):
            mapped_values = data_exp
        elif isinstance(data_exp, KM.Expression.ElementNonHistoricalExpression):
            neighbour_elems = KM.Expression.NodalNonHistoricalExpression(data_exp.GetModelPart())
            KOA.ContainerExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_elems)
            KOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, data_exp, neighbour_elems)
        elif isinstance(data_exp, KM.Expression.ConditionNonHistoricalExpression):
            neighbour_conds = KM.Expression.NodalNonHistoricalExpression(data_exp.GetModelPart())
            KOA.ContainerExpressionUtils.ComputeNumberOfNeighbourConditions(neighbour_conds)
            KOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, data_exp, neighbour_conds)

        filter_type = self._GetSolver().settings["filter_type"].GetString()
        if filter_type == "bulk_surface_shape" or filter_type == "general_vector":
            mapped_values.Evaluate(KOA.HELMHOLTZ_VECTOR_SOURCE)
        else:
            mapped_values.Evaluate(KOA.HELMHOLTZ_SCALAR_SOURCE)

    def _AssignNodalSolutionToDataExpression(self, output_data_exp_type):
        nodal_solution_field = KM.Expression.HistoricalExpression(self._GetComputingModelPart())
        filter_type = self._GetSolver().settings["filter_type"].GetString()
        if filter_type == "bulk_surface_shape" or filter_type == "general_vector":
            nodal_solution_field.Read(KOA.HELMHOLTZ_VECTOR)
        else:
            nodal_solution_field.Read(KOA.HELMHOLTZ_SCALAR)

        if output_data_exp_type == KM.Expression.NodalNonHistoricalExpression:
            non_hist_nodal_solution_field = KM.Expression.NodalNonHistoricalExpression(nodal_solution_field)
            return non_hist_nodal_solution_field
        elif output_data_exp_type == KM.Expression.ElementNonHistoricalExpression:
            mapped_elemental_solution_field = KM.Expression.ElementNonHistoricalExpression(self._GetComputingModelPart())
            KOA.ContainerExpressionUtils.MapNodalVariableToContainerVariable(mapped_elemental_solution_field, nodal_solution_field)
            return mapped_elemental_solution_field
        elif output_data_exp_type == KM.Expression.ConditionNonHistoricalExpression:
            mapped_condition_solution_field = KM.Expression.ConditionNonHistoricalExpression(self._GetComputingModelPart())
            KOA.ContainerExpressionUtils.MapNodalVariableToContainerVariable(mapped_condition_solution_field, nodal_solution_field)
            return mapped_condition_solution_field

    #### Public user interface functions ####
    def Initialize(self):
        if not self.initialized:
            super().Initialize()
            self._SetSolverMode()
            self.SetFilterRadius(self._GetSolver().settings["filter_radius"].GetDouble())
            self.SetBulkFilterRadius()
            self._SetHelmHoltzSourceMode()
            self.initialized = True

    def SetFilterRadius(self, filter_radius: float):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_RADIUS,filter_radius)

    def SetBulkFilterRadius(self):
        if self._GetSolver().settings["filter_type"].GetString() == "bulk_surface_shape":
            KOA.ImplicitFilterUtils.SetBulkRadiusForShapeFiltering(self._GetComputingModelPart())
        else:
            pass

    def RunSolver(self):
        self.Initialize()
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()

    def FilterField(self, unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        self._AssignDataExpressionToNodalSource(unfiltered_field)
        self.Initialize()
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()
        return self._AssignNodalSolutionToDataExpression(type(unfiltered_field))

    def FilterIntegratedField(self, unfiltered_field: KM.Expression.NodalNonHistoricalExpression) -> KM.Expression.HistoricalExpression:

        self._AssignDataExpressionToNodalSource(unfiltered_field)
        self.Initialize()
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(True)
        self.RunSolver()
        return self._AssignNodalSolutionToDataExpression(type(unfiltered_field))

    def UnFilterField(self, filtered_field: KM.Expression.NodalNonHistoricalExpression) -> KM.Expression.HistoricalExpression:

        self._AssignDataExpressionToNodalSource(filtered_field)
        self.Initialize()
        self._SetSolverMode(True)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()
        return self._AssignNodalSolutionToDataExpression(type(filtered_field))