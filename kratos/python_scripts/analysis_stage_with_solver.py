# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage

class AnalysisStageWithSolver(AnalysisStage):
    """
    Subclass of AnalysisStage that contains a solver.
    """
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

    def InitializeSolver(self):
        self._GetSolver().Initialize()
        #ModifyAfterSolverInitialize should be called here if needed

    def Finalize(self):
        super().Finalize()
        self._GetSolver().Finalize()

    def AddVariables(self):
        self._GetSolver().AddVariables()

    def ImportModelPart(self):
        self._GetSolver().ImportModelPart()

    def PrepareModelPart(self):
        self._GetSolver().PrepareModelPart()

    def AddDofs(self):
        self._GetSolver().AddDofs()

    def Predict(self):
        self._GetSolver().Predict()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self._GetSolver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()
        super().FinalizeSolutionStep()

    def SolveSolutionStep(self):
        return self._GetSolver().SolveSolutionStep()

    def AdvanceInTime(self):
        return self._GetSolver().AdvanceInTime(self.time)

    def Check(self):
        super().Check()
        self._GetSolver().Check()

    def Clear(self):
        super().Clear()
        self._GetSolver().Clear()

    def GetIsRestarted(self):
        return self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    def GetTime(self):
        return self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

    def SetTime(self, time):
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = time

    def GetStep(self):
        return self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP]

    def GetComputingModelPart(self):
        return self._GetSolver().GetComputingModelPart()

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")
