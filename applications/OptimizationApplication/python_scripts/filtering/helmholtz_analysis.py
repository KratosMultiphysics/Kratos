# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication.filtering.python_solvers_wrapper_implicit_filters as implicit_filter_solvers
import KratosMultiphysics.OptimizationApplication as KOA

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class HelmholtzAnalysis(AnalysisStage):
    """
    This class is the main-script of the implicit filtering analysis
    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases fail with a meaningful error message
        solver_settings = project_parameters["solver_settings"]
        if not solver_settings.Has("time_stepping"):
            raise Exception("Using the old way to pass the time_step, this was removed!")

        if not solver_settings.Has("domain_size"):
            raise Exception("Using the old way to pass the domain_size, this was removed!")

        if not solver_settings.Has("model_part_name"):
            raise Exception("Using the old way to pass the model_part_name, this was removed!")

        if not solver_settings.Has("echo_level"): # this is done to remain backwards-compatible
            raise Exception('"solver_settings" does not have "echo_level", please add it!')

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

    def _SetFilterRadius(self, filter_radius: float):
        self._GetComputingModelPart().GetProperties()[0].SetValue(KOA.HELMHOLTZ_RADIUS,filter_radius)

    def _SetHelmHoltzSourceMode(self, integrated_field=False):
        self._GetComputingModelPart().ProcessInfo.SetValue(KOA.HELMHOLTZ_INTEGRATED_FIELD, integrated_field)

    #### Public user interface functions ####
    def Initialize(self):
        if not self.initialized:
            super().Initialize()
            self.initialized = True

    def SetFilterRadius(self, filter_radius: float):
        self._GetComputingModelPart().GetProperties()[0].SetValue(KOA.HELMHOLTZ_RADIUS,filter_radius)

    def RunFilter(self):
        self.Initialize()
        self.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()

    def FilterField(self, unfiltered_field: KM.ContainerExpression.NodalNonHistoricalExpression) -> KM.ContainerExpression.HistoricalExpression:

        self.Initialize()

        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(False)

        unfiltered_field.Evaluate(KOA.HELMHOLTZ_SCALAR_SOURCE)

        self.RunFilter()

        filtered_field = KM.ContainerExpression.HistoricalExpression(self._GetComputingModelPart())
        filtered_field.Read(KOA.HELMHOLTZ_SCALAR)
        return filtered_field

    def FilterIntegratedField(self, unfiltered_field: KM.ContainerExpression.NodalNonHistoricalExpression) -> KM.ContainerExpression.HistoricalExpression:

        self.Initialize()

        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(True)

        unfiltered_field.Evaluate(KOA.HELMHOLTZ_SCALAR_SOURCE)
        self.RunFilter()

        filtered_field = KM.ContainerExpression.HistoricalExpression(self._GetSolver().GetComputingModelPart())
        filtered_field.Read(KOA.HELMHOLTZ_SCALAR)
        return filtered_field

    def UnFilterField(self, unfiltered_field: KM.ContainerExpression.NodalNonHistoricalExpression) -> KM.ContainerExpression.HistoricalExpression:

        self.Initialize()

        self._SetSolverMode(True)
        self._SetHelmHoltzSourceMode(False)

        unfiltered_field.Evaluate(KOA.HELMHOLTZ_SCALAR_SOURCE)
        self.RunFilter()

        filtered_field = KM.ContainerExpression.HistoricalExpression(self._GetSolver().GetComputingModelPart())
        filtered_field.Read(KOA.HELMHOLTZ_SCALAR)
        return filtered_field