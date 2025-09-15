# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class HighCycleFatigueAnalysis(StructuralMechanicsAnalysis):
    """This class is used to complement the structurea_mechanics_analysis
    when using the HCF constitutive law
    """

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            process = CLA.AdvanceInTimeHighCycleFatigueProcess(self._GetSolver().GetComputingModelPart(), self.project_parameters)
            process.Execute()
            time_incr = self._GetSolver().GetComputingModelPart().ProcessInfo[CLA.TIME_INCREMENT]
            self.time += time_incr
            self._GetSolver().GetComputingModelPart().ProcessInfo[CLA.TIME_INCREMENT] = 0.0
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()