from importlib import import_module
from KratosMultiphysics import python_solver
import KratosMultiphysics.NeuralNetworkApplication.input_dataclasses as InputDataclasses
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities
import KratosMultiphysics
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListDataWithLookback, ListNeuralNetworkData, NeuralNetworkData
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process_factory import NeuralNetworkProcessFactory

import tensorflow.keras as keras
import numpy as np

def CreateSolver(project_parameters, model = None):
    return NeuralNetworkSolver(project_parameters,model)

class NeuralNetworkSolver(PythonSolver):
    """
    This class is the based on the AnalysisStage from Kratos but regarding the creation and training of the neural
    network model.
    """

    # TODO: Addapt solver to PythonSolver structure (project parameters)
    def __init__(self, project_parameters, model = None):      

        #TODO: This should be done through default parameters

        self.echo_level = project_parameters["echo_level"].GetInt()

        self.model = model
        self.model_geometry_name = project_parameters["model_part_name"].GetString()
        self.model_geometry_file = project_parameters["model_part_file"].GetString()
        self.model_import_settings = KratosMultiphysics.Parameters()
        self.model_import_settings.AddEmptyValue("input_type")
        self.model_import_settings["input_type"].SetString("mdpa")
        self.model_import_settings.AddEmptyValue("input_filename")
        self.model_import_settings["input_filename"].SetString(self.model_geometry_file)
        self.time_buffer = project_parameters["time_buffer"].GetInt()
        self.timestep = project_parameters["timestep"].GetDouble()
        self.end_time = project_parameters["end_time"].GetDouble()
        self.start_time = project_parameters["start_time"].GetDouble()
        try:
            self.record = project_parameters["record"].GetBool()
        except RuntimeError:
            self.record = False
        try:
            self.external_model = project_parameters["external_model"].GetBool()
        except RuntimeError:
            self.external_model = False
        try:
            self.only_input = project_parameters["only_input"].GetBool()
        except RuntimeError:
            self.only_input = False
        try:
            self.timesteps_as_features = project_parameters["timesteps_as_features"].GetBool()
        except RuntimeError:
            self.timesteps_as_features = False
        try:
            self.feaures_as_timestpes = project_parameters["features_as_timesteps"].GetBool()
        except RuntimeError:
            self.feaures_as_timestpes = False
        try:
            self.soft_start_flag = project_parameters["soft_start"].GetBool()
        except RuntimeError:
            self.soft_start_flag = True
        try:
            self.dimension_in = project_parameters["dimension_input"].GetInt()
        except RuntimeError:
            self.dimension_in = 1
        try:
            self.reorder_partitions = project_parameters["reorder_partitions"].GetInt()
        except RuntimeError:
            self.reorder_partitions = 1

        output_var_names = project_parameters["output_variables"]
        output_sources_names = project_parameters["output_sources"]
        output_variable_names = [ output_var_names[i].GetString() for i in range( output_var_names.size() ) ]
        self.output_sources = [ output_sources_names[i].GetString() for i in range( output_sources_names.size() ) ]
        self.output_variables = [ KratosMultiphysics.KratosGlobals.GetVariable( var ) for var in output_variable_names ]

    def Initialize(self):
        self.LoadGeometry()
        self.input_data_structure = InputDataclasses.NeuralNetworkData()

        self.preprocessed_data_structure = InputDataclasses.NeuralNetworkData()
        self.preprocessed_previous = InputDataclasses.NeuralNetworkData()
        
        self.output_data_structure = InputDataclasses.NeuralNetworkData()
        self.time = 0.0


    def InitializeSolutionStep(self):
        model_input_value_list = []
        for node in self.model_geometry.Nodes:
            model_input_value_list.append([node.X, node.Y, self.time])
    
        self.input_from_modelpart = np.array(model_input_value_list)

    def SolveSolutionStep(self):

        print("Receiving data into the neural network model")

        # Retrieve input from geometry
        self.input_data_structure.UpdateData(self.input_from_modelpart)

        # Predict (and invert the transformations) from the new input and update it to the output
        self.output_data_structure.UpdateData(self.PredictNeuralNetwork(data_structure_in = self.input_data_structure))
        output_value_list = np.squeeze(self.output_data_structure.ExportAsArray())

        if self.time >= self.time_buffer: # TODO: Check for consistency with timesteps nad time
            output_value_index = 0

            for i, variable in enumerate(self.output_variables):
                for node, node_id in zip(self.model_geometry.Nodes, range(self.model_geometry.NumberOfNodes())):
                    node.SetSolutionStepValue(variable, 0, output_value_list[node_id,i])
                    output_value_index += node_id


    def PredictNeuralNetwork(self, data_structure_in = None):
       
        data_out = NeuralNetworkData()
        self._PrintInfo("Predicting with the Physics Informed Neural Network...")

        if self.external_model:
            output = self.neural_network_model.predict(data_structure_in.ExportAsArray())
            if not output is None:
                data_out.UpdateData(output)
        else:
            for process in self._GetListOfProcesses():
                try:
                    output = process.Predict(self.neural_network_model, data_structure_in)
                    if not output is None:
                        data_out.UpdateData(output)
                except AttributeError:
                    pass
        return data_out

    def FinalizeSolutionStep(self):
        pass 
    
    def AdvanceInTime(self, previous_time):
        self.time = previous_time + self.timestep
        if hasattr(self, 'model_geometry'):
            new_time = self.model_geometry.ProcessInfo[KratosMultiphysics.TIME] + self.timestep
            self.model_geometry.ProcessInfo.SetValue(KratosMultiphysics.TIME, new_time)
            new_step = self.model_geometry.ProcessInfo[KratosMultiphysics.STEP] + 1
            self.model_geometry.ProcessInfo.SetValue(KratosMultiphysics.STEP, new_step)
        return self.time

    def LoadNeuralNetworkModel(self, neural_network_model):
        self.neural_network_model = neural_network_model
    
    def LoadGeometry(self):
        self.model_geometry = self.model[self.model_geometry_name]

    def ImportModelPart(self):
        self._ImportModelPart(self.model[self.model_geometry_name], self.model_import_settings)

    def LoadProcessesList(self,list_of_processes):
        self._list_of_processes = list_of_processes

    def _GetListOfProcesses(self):
        return self._list_of_processes
    
    def _PrintInfo(self, message):
        """This function prints info messages if the echo level is greater than 1.
        """
        if self.echo_level > 0:
            print(message)   