# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication.filtering.python_solvers_wrapper_implicit_filters as implicit_filter_solvers
from KratosMultiphysics.OptimizationApplication.filtering.helmholtz_bulk_surface_solver import HelmholtzBulkSurfaceSolver
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
import KratosMultiphysics.OptimizationApplication as KOA

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class HelmholtzAnalysis(AnalysisStage):
    """
    This class is the main-script of the implicit filtering analysis
    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):

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

    def _SetSolverMode(self, invers_mode=False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.COMPUTE_HELMHOLTZ_INVERSE, invers_mode)

    def _SetHelmHoltzSourceMode(self, integrated_field=False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_INTEGRATED_FIELD, integrated_field)

    def _MapDataExpressionToNodalExpression(self, data_exp: ContainerExpressionTypes):
        if type(data_exp) == KM.ContainerExpression.NodalNonHistoricalExpression:
            return data_exp
        elif type(data_exp) == KM.ContainerExpression.ElementNonHistoricalExpression:
            neighbour_elems = KM.ContainerExpression.NodalNonHistoricalExpression(data_exp.GetModelPart())
            KOA.ContainerExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_elems)
            mapped_values = KM.ContainerExpression.NodalNonHistoricalExpression(data_exp.GetModelPart())
            KOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, data_exp, neighbour_elems)
            return mapped_values
        elif type(data_exp) == KM.ContainerExpression.ConditionNonHistoricalExpression:
            neighbour_conds = KM.ContainerExpression.NodalNonHistoricalExpression(data_exp.GetModelPart())
            KOA.ContainerExpressionUtils.ComputeNumberOfNeighbourConditions(neighbour_conds)
            mapped_values = KM.ContainerExpression.NodalNonHistoricalExpression(data_exp.GetModelPart())
            KOA.ContainerExpressionUtils.MapContainerVariableToNodalVariable(mapped_values, data_exp, neighbour_conds)
            return mapped_values

    def _MapNodalSolutionToDataExpression(self, output_data_exp_type):

        nodal_solution_field = KM.ContainerExpression.HistoricalExpression(self._GetComputingModelPart())
        if type(self._GetSolver()) == HelmholtzBulkSurfaceSolver:
            nodal_solution_field.Read(KOA.HELMHOLTZ_VECTOR)
        else:
            nodal_solution_field.Read(KOA.HELMHOLTZ_SCALAR)

        if output_data_exp_type == KM.ContainerExpression.NodalNonHistoricalExpression:
            non_hist_nodal_solution_field = KM.ContainerExpression.NodalNonHistoricalExpression(nodal_solution_field)
            return non_hist_nodal_solution_field
        elif output_data_exp_type == KM.ContainerExpression.ElementNonHistoricalExpression:
            mapped_elemental_solution_field = KM.ContainerExpression.ElementNonHistoricalExpression(self._GetComputingModelPart())
            KOA.ContainerExpressionUtils.MapNodalVariableToContainerVariable(mapped_elemental_solution_field, nodal_solution_field)
            return mapped_elemental_solution_field
        elif output_data_exp_type == KM.ContainerExpression.ConditionNonHistoricalExpression:
            mapped_condition_solution_field = KM.ContainerExpression.ConditionNonHistoricalExpression(self._GetComputingModelPart())
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
        if type(self._GetSolver()) == HelmholtzBulkSurfaceSolver:
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

        mapped_unfiltered_field = self._MapDataExpressionToNodalExpression(unfiltered_field)

        size = KOA.ContainerExpressionUtils.GetFlattenedSize(mapped_unfiltered_field)
        if size==3 and type(self._GetSolver()) == HelmholtzBulkSurfaceSolver:
            self.BulkSurfaceFilterVectorField(mapped_unfiltered_field)
        elif size==3:
            self.FilterVectorField(mapped_unfiltered_field)
        elif size==1:
            self.FilterScalarField(mapped_unfiltered_field)

        return self._MapNodalSolutionToDataExpression(type(unfiltered_field))

    def FilterScalarField(self, mapped_unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        mapped_unfiltered_field.Evaluate(KOA.HELMHOLTZ_SCALAR_SOURCE)
        self.Initialize()
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()

    def BulkSurfaceFilterVectorField(self, mapped_unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        mapped_unfiltered_field.Evaluate(KOA.HELMHOLTZ_VECTOR_SOURCE)
        self.Initialize()
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()

    def FilterVectorField(self, mapped_unfiltered_field: ContainerExpressionTypes) -> ContainerExpressionTypes:

        filtered_field = KM.ContainerExpression.NodalNonHistoricalExpression(self._GetComputingModelPart())
        filtered_field.SetData(KM.Array3([0, 0, 0]))

        for i in range(3):
            KOA.ImplicitFilterUtils.AssignVectorNodalExpressionToScalarVariable(KOA.HELMHOLTZ_SCALAR_SOURCE,mapped_unfiltered_field,i)
            self.Initialize()
            self._SetSolverMode(False)
            self._SetHelmHoltzSourceMode(False)
            self.RunSolver()
            KOA.ImplicitFilterUtils.AssignScalarVariableToVectorNodalExpression(KOA.HELMHOLTZ_SCALAR_SOURCE,filtered_field,i)

    def FilterIntegratedField(self, unfiltered_field: KM.ContainerExpression.NodalNonHistoricalExpression) -> KM.ContainerExpression.HistoricalExpression:

        self._MapDataExpressionToNodalSource(unfiltered_field)

        self.Initialize()
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(True)
        self.RunSolver()

        return self._MapNodalSolutionToDataExpression(type(unfiltered_field))



    def UnFilterField(self, filtered_field: KM.ContainerExpression.NodalNonHistoricalExpression) -> KM.ContainerExpression.HistoricalExpression:

        self._MapDataExpressionToNodalSource(filtered_field)

        self.Initialize()
        self._SetSolverMode(True)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()

        return self._MapNodalSolutionToDataExpression(type(filtered_field))