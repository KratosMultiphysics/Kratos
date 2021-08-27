import KratosMultiphysics as KM

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing NeuralNetworkApplication
if not CheckIfApplicationsAvailable("NeuralNetworkApplication"):
    raise ImportError("The NeuralNetworkApplication is not available!")
from KratosMultiphysics.NeuralNetworkApplication.neural_network_analysis import NeuralNetworkAnalysis

def Create(settings, model, solver_name):
    return NeuralNetworkWrapper(settings, model, solver_name)

class NeuralNetworkWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the NeuralNetworkApplication of Kratos"""
    def __init__(self, settings, model, solver_name):
        super(NeuralNetworkWrapper, self).__init__(settings, model, solver_name)

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        self.input_variable = self.settings["solver_wrapper_settings"]["input_variable"].GetString()
        self.output_variable = self.settings["solver_wrapper_settings"]["output_variable"].GetString()
        self.mp = self.model.CreateModelPart("NeuralNetwork")
        # self.mp.AddNodalSolutionStepVariable(getattr(KM, self.input_variable))
        self.mp.AddNodalSolutionStepVariable(getattr(KM, "FORCE_Y"))
        # self.mp.AddNodalSolutionStepVariable(getattr(KM, self.output_variable))
        self.mp.AddNodalSolutionStepVariable(getattr(KM, "DISPLACEMENT_Y"))
        mdpa_file_name = self.settings["solver_wrapper_settings"]["mdpa_file_name"].GetString()
    
        # KM.ModelPartIO(mdpa_file_name).ReadModelPart(self.mp)

        # self.mp.ProcessInfo[KM.DOMAIN_SIZE] = self.project_parameters["domain_size"].GetInt()

    def AdvanceInTime(self, current_time):
        return 0.0

    def SolveSolutionStep(self):
        with self.thread_manager:
            print("Receiving data into the neural network model")
            input_value =[]
            # for variable in self.dict_input.items():
            #         input_value.append(self.GetInterfaceData(variable).GetData())

            output_value =[]
            # for variable in self.dict_output.items():
            #         output_value.append(self.GetInterfaceData(variable).SetData())

            input_value.append(self.GetInterfaceData(self.input_variable).GetData())
            self._analysis_stage.Initialize(data_input = input_value[0])
            output_value.append(self._analysis_stage.Predict())
            self.GetInterfaceData(self.output_variable).SetData(output_value[0])
            print("Predicted data using neural network model.")

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def OutputSolutionStep(self):
        pass

    def _CreateAnalysisStage(self):
        return NeuralNetworkAnalysis(self.project_parameters)