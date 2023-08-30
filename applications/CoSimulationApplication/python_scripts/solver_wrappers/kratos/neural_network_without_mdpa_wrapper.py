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
        # TODO: The input and output can only be one variable in the current state.
        self.input_variable = self.settings["solver_wrapper_settings"]["input_variable"].GetString()
        self.output_variable = self.settings["solver_wrapper_settings"]["output_variable"].GetString()
        self.mp = self.model.CreateModelPart("NeuralNetwork")
        self.lookback = self.settings["solver_wrapper_settings"] ["lookback"].GetInt()
        try:
            self.record = settings["solver_wrapper_settings"]["record"].GetBool()
        except RuntimeError:
            self.record = False
        try:
            self.timesteps_as_features = settings["solver_wrapper_settings"]["timesteps_as_features"].GetBool()
        except RuntimeError:
            self.timesteps_as_features = False
        try:
            self.feaures_as_timestpes = settings["solver_wrapper_settings"]["features_as_timesteps"].GetBool()
        except RuntimeError:
            self.feaures_as_timestpes = False
        try:
            self.soft_start_flag = settings["solver_wrapper_settings"]["soft_start"].GetBool()
        except RuntimeError:
            self.soft_start_flag = True

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
                if self.record:
                    self.preprocessed_data_structure.CheckRecordAndUpdate(self.output_data_structure.ExportAsArray())
            # Preprocess input and output 
            [self.input_data_structure, self.output_data_structure] = self._analysis_stage.Preprocessing(
                data_in = self.input_data_structure, data_out = self.output_data_structure)
            # Update input to the structure with the new preprocessed input
            if self.record:
                self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                # This flag softens the initial lookback in the record (stabilizes the behaviour)
                if self.soft_start_flag:
                    for i in range(self.lookback-1):
                        self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                    self.soft_start_flag = False
                self.record_data = self.preprocessed_data_structure.data
            else:
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
        if self.record:
            self.preprocessed_data_structure.data=self.record_data
        self.preprocessed_previous =self.preprocessed_data_structure

    def OutputSolutionStep(self):
        pass

    def _CreateAnalysisStage(self):
        return NeuralNetworkAnalysis(self.project_parameters)