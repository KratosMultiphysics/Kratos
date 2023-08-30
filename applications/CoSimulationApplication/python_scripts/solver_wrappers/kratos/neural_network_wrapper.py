# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing NeuralNetwork
if not CheckIfApplicationsAvailable("NeuralNetworkApplication"):
    raise ImportError("The NeuralNetworkApplication is not available!")
from KratosMultiphysics.NeuralNetworkApplication.neural_network_analysis import NeuralNetworkAnalysis

def Create(settings, model, solver_name):
    return NeuralNetworkWrapper(settings, model, solver_name)

class NeuralNetworkWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the NeuralNetworkApplication of Kratos"""
    def _CreateAnalysisStage(self):
        return NeuralNetworkAnalysis(self.project_parameters, self.model)
   
    def SolveSolutionStep(self):
        with self.thread_manager:
            self._analysis_stage._GetNeuralNetworkSolver().SolveSolutionStep()
    
    def Predict(self):
        with self.thread_manager:
            self._analysis_stage._GetNeuralNetworkSolver().Predict()
    
    def AdvanceInTime(self, current_time):
        with self.thread_manager:
            new_time = self._analysis_stage._GetNeuralNetworkSolver().AdvanceInTime(current_time)
        self._analysis_stage.time = new_time # only needed to print the time correctly
        return new_time