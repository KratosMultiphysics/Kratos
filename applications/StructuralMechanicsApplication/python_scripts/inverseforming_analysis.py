# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class InverseFormingAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)
    
    def ModifyInitialGeometry(self):
        KratosMultiphysics.Logger.PrintInfo("::[InverseFormingAnalysis]:: ", "Initializing guess with vertical projection")
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            node.Z0 = 0.0
        KratosMultiphysics.Logger.PrintInfo("::[InverseFormingAnalysis]:: ", "InitialGeometry flattened")
         # *CHECK* Correct place to calculate normals? Possible in element, strategy or else?
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()
        normal_calculation_utils.CalculateUnitNormalsNonHistorical(self._GetSolver().GetComputingModelPart(), 0)
        KratosMultiphysics.Logger.PrintInfo("::[InverseFormingAnalysis]:: ", "UnitNormals Calculated")
    
    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        mysolver = self._GetSolver()
        print(mysolver.settings)
        self.originalZ = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            self.originalZ.append(node.Z)
        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            self.ScalingGeometry()
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
    
    def ScalingGeometry(self):
        self.scaling = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME]
        current_step = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP]
        print("scaling:", self.scaling)
        print("current step:", current_step)
        print("step scaling:", self.scaling * current_step)
        print("current time: ", self.time)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            node.Z = self.originalZ[node.Id - 1] * self.scaling * current_step
    
    def CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        # Wrapper to call the mangled method
        return self._AnalysisStage__CheckIfSolveSolutionStepReturnsAValue(is_converged)