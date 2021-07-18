import KratosMultiphysics as KM

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing StructuralMechanics
if not CheckIfApplicationsAvailable("NeuralNetworkApplication"):
    raise ImportError("The NeuralNetworkApplication is not available!")
from KratosMultiphysics.NeuralNetworkApplication.neural_network_analysis import NeuralNetworkAnalysis

def Create(settings, solver_name):
    return NeuralNetworkWrapper(settings, solver_name)

class NeuralNetworkWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the NeuralNetworkApplication of Kratos"""
    def __init__(self, settings, solver_name):
        super(NeuralNetworkWrapper, self).__init__(settings, solver_name)

    def AdvanceInTime(self, current_time):
        return 0.0

    def SolveSolutionStep(self):

        print("Receiving data into the neural network model")
        input_value =[]
        for variable in self.dict_input.items():
                input_value.append(self.GetInterfaceData(variable).GetData())

        output_value =[]
        for variable in self.dict_output.items():
                output_value.append(self.SetInterfaceData(variable).GetData())

        print("Predicted data using neural network model")

    def _CreateAnalysisStage(self):
        return NeuralNetworkAnalysis(self.project_parameters)