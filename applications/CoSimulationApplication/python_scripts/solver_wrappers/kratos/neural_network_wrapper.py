import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication.input_dataclasses as InputDataclasses

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
        self.lookback = self.settings["solver_wrapper_settings"] ["lookback"].GetInt()
    
        # KM.ModelPartIO(mdpa_file_name).ReadModelPart(self.mp)

        # self.mp.ProcessInfo[KM.DOMAIN_SIZE] = self.project_parameters["domain_size"].GetInt()

    def AdvanceInTime(self, current_time):
        return 0.0

    def SolveSolutionStep(self):
        with self.thread_manager:
            print("Receiving data into the neural network model")

            # Restarting preprocessed data if no FinalizeSolutionStep took place
            try:
                self.preprocessed_data_structure = self.preprocessed_previous
            except AttributeError:
                pass

            # Retrieve input from interface
            self.input_data_structure.UpdateData(self.GetInterfaceData(self.input_variable).GetData())
            # Initialize output from interface in first iteration
            if self.output_data_structure.data == None:
                self.output_data_structure.UpdateData(self.GetInterfaceData(self.output_variable).GetData())
            # Initialize preprocessed input to the network in first iteration
            if self.lookback>0 and not self.preprocessed_data_structure.lookback_state:
                self.preprocessed_data_structure.CheckLookbackAndUpdate(self.output_data_structure.ExportAsArray())
            # Preprocess input and output 
            [self.input_data_structure, self.output_data_structure] = self._analysis_stage.Preprocessing(
                data_in = self.input_data_structure, data_out = self.output_data_structure)
            # Update input to the structure with the new preprocessed input
            self.preprocessed_data_structure.UpdateData(self.input_data_structure.ExportAsArray())
            # Predict (and invert the transformations) from the new input and update it to the output
            self.output_data_structure.UpdateData(self._analysis_stage.Predict(data_structure_in = self.preprocessed_data_structure))
            # Export output to the interface
            self.GetInterfaceData(self.output_variable).SetData(self.output_data_structure.ExportAsArray())
            print("Predicted data using neural network model.")
            

    def Initialize(self):
        self.input_data_structure = InputDataclasses.NeuralNetworkData()
        if self.lookback>0:
            self.preprocessed_data_structure = InputDataclasses.DataWithLookback(lookback_index=self.lookback)  
        else:
            self.preprocessed_data_structure = InputDataclasses.NeuralNetworkData()
        self.output_data_structure = InputDataclasses.NeuralNetworkData()

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        
        pass

    def Predict(self):
        pass

    def FinalizeSolutionStep(self):
        self.preprocessed_previous =self.preprocessed_data_structure
        pass

    def OutputSolutionStep(self):
        pass

    def _CreateAnalysisStage(self):
        return NeuralNetworkAnalysis(self.project_parameters)