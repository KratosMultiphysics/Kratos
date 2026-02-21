from typing import Any, Optional

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
import KratosMultiphysics.OptimizationApplication.filtering.python_solvers_wrapper_implicit_filters as implicit_filter_solvers

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

class HelmholtzAnalysis(AnalysisStage):
    """
    This class is the main-script of the implicit filtering analysis
    It can be imported and used as "black-box"
    """
    def __init__(self, model: KM.Model, project_parameters: KM.Parameters):
        super().__init__(model, project_parameters)
        self.__source_data: 'Optional[KM.TensorAdaptors.DoubleTensorAdaptor]' = None
        self.__neighbour_entities: 'dict[Any, KM.TensorAdaptors.IntTensorAdaptor]' = {}

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

    def FilterField(self, unfiltered_field: KM.TensorAdaptors.DoubleTensorAdaptor) -> KM.TensorAdaptors.DoubleTensorAdaptor:

        self.__AssignTensorDataToNodalSource(unfiltered_field)
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()
        return self.__AssignNodalSolutionToTensorData()

    def FilterIntegratedField(self, unfiltered_field: KM.TensorAdaptors.DoubleTensorAdaptor) -> KM.TensorAdaptors.DoubleTensorAdaptor:

        self.__AssignTensorDataToNodalSource(unfiltered_field)
        self._SetSolverMode(False)
        self._SetHelmHoltzSourceMode(True)
        self.RunSolver()
        return self.__AssignNodalSolutionToTensorData()

    def UnFilterField(self, filtered_field: KM.TensorAdaptors.DoubleTensorAdaptor) -> KM.TensorAdaptors.DoubleTensorAdaptor:

        self.__AssignTensorDataToNodalSource(filtered_field)
        self._SetSolverMode(True)
        self._SetHelmHoltzSourceMode(False)
        self.RunSolver()
        return self.__AssignNodalSolutionToTensorData()

    def AssignTensorDataToNodalSolution(self, data_exp: KM.TensorAdaptors.DoubleTensorAdaptor) -> None:
        ta = KM.TensorAdaptors.DoubleTensorAdaptor(data_exp, copy=False)
        if not isinstance(ta.GetContainer(), KM.NodesArray):
            neighbours = self.__GetNeighbourEntities(data_exp)
            ta = KOA.OptimizationUtils.MapContainerDataToNodalData(data_exp, neighbours.GetContainer())
        KM.TensorAdaptors.HistoricalVariableTensorAdaptor(ta, self._GetSolver().GetSolvingVariable(), copy=False).StoreData()

    def __AssignTensorDataToNodalSource(self, data_exp: KM.TensorAdaptors.DoubleTensorAdaptor):
        self.__source_data = KM.TensorAdaptors.DoubleTensorAdaptor(data_exp)

        # it is better to work on the model part of the data_exp rather than the internal
        # model part created with ConnectivityPreserveModelPart because, then all the outputs will
        # be written to one vtu file automatically, otherwise, there will be outputs distributed among
        # the helmholtz model part and the original model part. This is safer because the helmholtz modelpart
        # is created using the ConnectivityPreserveModeller which preserves the same nodes and condition/element
        # data containers.

        ta = KM.TensorAdaptors.DoubleTensorAdaptor(data_exp, copy=False)
        if not isinstance(ta.GetContainer(), KM.NodesArray):
            neighbours = self.__GetNeighbourEntities(data_exp)
            ta = KOA.OptimizationUtils.MapContainerDataToNodalData(data_exp, neighbours.GetContainer())
        KM.TensorAdaptors.VariableTensorAdaptor(ta, self._GetSolver().GetSourceVariable(), copy=False).StoreData()

    def __AssignNodalSolutionToTensorData(self) -> KM.TensorAdaptors.DoubleTensorAdaptor:
        if self.__source_data is None:
            raise RuntimeError("The __AssignTensorDataToNodalSource should be called first.")

        # it is better to work on the model part of the data_exp rather than the internal
        # model part created with ConnectivityPreserveModelPart because, then all the outputs will
        # be written to one vtu file automatically, otherwise, there will be outputs distributed among
        # the helmholtz model part and the original model part. This is safer because the helmholtz modelpart
        # is created using the ConnectivityPreserveModeller which preserves the same nodes and condition/element
        # data containers.
        neighbours = self.__GetNeighbourEntities(self.__source_data)
        solution = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(neighbours.GetContainer(), self._GetSolver().GetSolvingVariable())
        solution.CollectData()

        if isinstance(self.__source_data.GetContainer(), KM.NodesArray):
            return KM.TensorAdaptors.DoubleTensorAdaptor(solution, copy=False)
        else:
            return KOA.OptimizationUtils.MapNodalDataToContainerData(solution, self.__source_data.GetContainer(), neighbours)

    def __GetNeighbourEntities(self, data_exp: KM.TensorAdaptors.DoubleTensorAdaptor) -> KM.TensorAdaptors.IntTensorAdaptor:
        # following makes the number of neighbours computation to be executed once
        # per given container, hence if the mesh element/connectivity changes
        # this computation needs to be redone. Especially in the case if MMG is
        # used for re-meshing.
        key = data_exp.GetContainer()
        if key not in  self.__neighbour_entities.keys():
            if not isinstance(key, KM.NodesArray):
                ta = KM.TensorAdaptors.NodalNeighbourCountTensorAdaptor(self._GetSolver().GetComputingModelPart().Nodes, data_exp.GetContainer())
                ta.CollectData()
                self.__neighbour_entities[key] = ta
            else:
                ta = KM.TensorAdaptors.IntTensorAdaptor(key, KM.IntNDData([len(key)]), copy=False)
                ta.data[:] = 1
                self.__neighbour_entities[key] = ta

        return self.__neighbour_entities[key]