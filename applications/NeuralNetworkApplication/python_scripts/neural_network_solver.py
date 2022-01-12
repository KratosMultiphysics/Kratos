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


        self.model = model
        self.model_geometry_name = project_parameters["model_part_name"].GetString()


        self.lookback = project_parameters ["lookback"].GetInt()
        self.time_buffer = project_parameters["time_buffer"].GetInt()
        self.timestep = project_parameters["timestep"].GetDouble()
        self.end_time = project_parameters["end_time"].GetDouble()
        self.start_time = project_parameters["start_time"].GetDouble()
        try:
            self.record = project_parameters["record"].GetBool()
        except RuntimeError:
            self.record = False
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

        # TODO: This block is shared with data_generator_process, it could be separated and shared throug a function

        # getting the ModelPart from the Model
        self.output_model_part_name = project_parameters["output_model_part"].GetString()
        if self.output_model_part_name == "":
            raise Exception('No "output_model_part" was specified!')

        # getting the input ModelPart from the Model
        self.input_model_part_name = project_parameters["input_model_part"].GetString()
        if self.input_model_part_name == "":
            raise Exception('No "input_model_part" was specified!')

        # retrieving the input variables
        input_var_names = project_parameters["input_variables"]
        input_sources_names = project_parameters["input_sources"]
        variable_names = [ input_var_names[i].GetString() for i in range( input_var_names.size() ) ]
        self.input_sources = [ input_sources_names[i].GetString() for i in range( input_sources_names.size() ) ]
        self.input_variables = [ KratosMultiphysics.KratosGlobals.GetVariable( var ) for var in variable_names ]
        if len(self.input_variables) == 0:
            raise Exception('No variables specified for input!')
        if not (len(self.input_variables) == len(self.input_sources)):
            raise Exception('The number of input variables and sources are different.')
        self.dict_input = dict(zip(self.input_variables, self.input_sources))
        # getting input order
        try:
            self.input_order = project_parameters["input_order"].GetString()
        except RuntimeError:
            self.input_order = "variables_first"

        # retrieving the output variables
        output_var_names = project_parameters["output_variables"]
        output_sources_names = project_parameters["output_sources"]
        variable_names = [ output_var_names[i].GetString() for i in range( output_var_names.size() ) ]
        self.output_sources = [ output_sources_names[i].GetString() for i in range( output_sources_names.size() ) ]
        self.output_variables = [ KratosMultiphysics.KratosGlobals.GetVariable( var ) for var in variable_names ]
        if len(self.output_variables) == 0:
            raise Exception('No variables specified for output!')
        if not (len(self.output_variables) == len(self.output_sources)):
            raise Exception('The number of output variables and sources are different.')
        self.dict_output = dict(zip(self.output_variables, self.output_sources))
        # getting output order
        try:
            self.output_order = project_parameters["output_order"].GetString()
        except RuntimeError:
            self.output_order = "variables_first"

    def Initialize(self):
        self.LoadGeometry()

        self.input_data_structure = InputDataclasses.NeuralNetworkData()
        if self.lookback>0:
            self.preprocessed_data_structure = InputDataclasses.DataWithLookback(lookback_index=self.lookback)  
            self.preprocessed_previous = InputDataclasses.DataWithLookback(lookback_index=self.lookback)
        else:
            self.preprocessed_data_structure = InputDataclasses.NeuralNetworkData()
            self.preprocessed_previous = InputDataclasses.NeuralNetworkData()
        self.output_data_structure = InputDataclasses.NeuralNetworkData()
        self.time = 0.0
        self.output_model_part = self.model[self.output_model_part_name]
        self.input_model_part = self.model[self.input_model_part_name]


    def InitializeSolutionStep(self):

        input_value_list=[]

        for variable, source in self.dict_input.items():
            # Process related variables (e.g. TIME, STEP)
            if source == 'process':
                input_value_list.append(self.input_model_part.ProcessInfo[variable])
            # Node properties (e.g. position)
            elif source == 'node':
                for node in self.input_model_part.Nodes:
                    input_value_list.append(getattr(node,variable.Name()))
            # Node step values (e.g. variables like displacement)
            elif source == "solution_step":
                for node in self.input_model_part.Nodes:
                    input_value_list.append(node.GetSolutionStepValue(variable,0))
            # Condition values
            elif source == "condition":
                for condition in self.input_model_part.GetConditions():
                    input_value_list.append(condition.GetValue(variable))
        # Reorder if indicated
        if self.input_order == 'sources_first':
            try:
                input_value_list = self._OrderSourcesFirst(input_value_list, self.dict_output.items())
            except IndexError:
                pass

        self.input_from_modelpart = np.array(input_value_list)


    def SolveSolutionStep(self):

        print("Receiving data into the neural network model")

        # Restarting preprocessed data if no FinalizeSolutionStep took place
        try:
            if self.converge_flag == False:
                self.preprocessed_data_structure.UpdateData(self.preprocessed_previous.data)
                if hasattr(self.preprocessed_previous,'lookback_data'):
                    if not self.only_input:
                        self.preprocessed_data_structure.UpdateLookbackAll(self.preprocessed_previous.lookback_data)
                elif self.preprocessed_data_structure.data is None:
                    if not self.only_input:
                        self.preprocessed_data_structure.UpdateLookbackAll(np.zeros_like(self.preprocessed_data_structure.lookback_data))
                    self.soft_start_flag = True
                    if self.only_input:
                        if self.record:
                            self.preprocessed_data_structure.record_data = False
                        self.preprocessed_data_structure.lookback_state = False
                        delattr(self.preprocessed_data_structure, 'lookback_data')
                else:
                    if self.only_input:
                        delattr(self.preprocessed_data_structure, 'lookback_data')
                        self.preprocessed_data_structure.lookback_state = False
                    else:
                        self.preprocessed_data_structure.lookback_state = False
                        if self.record:
                            self.preprocessed_data_structure.record_data = False
            else: self.converge_flag = False
        except AttributeError:
            pass

        # Retrieve input from geometry
        self.input_from_modelpart = np.squeeze(np.reshape(self.input_from_modelpart, (int(self.input_from_modelpart.size/self.dimension_in), self.dimension_in)))
        self.input_data_structure.UpdateData(self.input_from_modelpart)

        # Initialize output from interface in first iteration
        if self.output_data_structure.data is None:
            for variable, source in self.dict_output.items():
                output_value_list = []
                # Process related variables (e.g. TIME, STEP)
                if source == 'process':
                    output_value_list.append(self.output_model_part.ProcessInfo[variable])
                # Node properties (e.g. position)
                elif source == 'node':
                    for node in self.output_model_part.Nodes:
                        output_value_list.append(getattr(node,variable.Name()))
                # Node step values (e.g. variables like displacement)
                elif source == "solution_step":
                    for node in self.output_model_part.Nodes:
                        output_value_list.append(node.GetSolutionStepValue(variable,0))
                # Condition values
                elif source == "condition":
                    for condition in self.output_model_part.GetConditions():
                        output_value_list.append(condition.GetValue(variable))
            # Reorder if indicated
            if self.output_order == 'sources_first':
                try:
                    output_value_list = self._OrderSourcesFirst(output_value_list, self.dict_output.items())
                except IndexError:
                    pass
            self.output_data_structure.UpdateData(output_value_list)
        
        # Preprocess input and output 
        [self.input_data_structure, self.output_data_structure] = self.Preprocessing(
            data_in = self.input_data_structure, data_out = self.output_data_structure)
        # Initialize preprocessed input to the network in first iteration
        if self.lookback>0 and not self.preprocessed_data_structure.lookback_state:
            if not self.only_input:
                self.preprocessed_data_structure.CheckLookbackAndUpdate(self.output_data_structure.ExportAsArray())
            if self.record:
                self.preprocessed_data_structure.CheckRecordAndUpdate(self.input_data_structure.ExportAsArray())
        # Update input to the structure with the new preprocessed input
        if self.record:
            if self.preprocessed_data_structure.data is None:
                self.preprocessed_data_structure.record_data = False
                self.preprocessed_data_structure.CheckRecordAndUpdate(self.input_data_structure.ExportAsArray())
            # This flag softens the initial lookback in the record (stabilizes the behaviour)
            if self.soft_start_flag:
                for i in range(self.lookback-1):
                    if not self.only_input:
                        self.preprocessed_data_structure.CheckLookbackAndUpdate(self.output_data_structure.ExportAsArray())
                    self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                self.soft_start_flag = False
            else:
                self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
            self.record_data = self.preprocessed_data_structure.data
        else:
            self.preprocessed_data_structure.UpdateData(self.input_data_structure.ExportAsArray())
        # Set the reordering of the preprocessed data structure
        if self.reorder_partitions > 1:
            self.preprocessed_data_structure.reorder_partitions = self.reorder_partitions
            
            # Undo the reorder if indicated
        if self.output_order == 'sources_first':
            try:
                output_value_list = self._OrderVariablesFirst(output_value_list, self.dict_output.items())
            except IndexError:
                pass

        # Predict (and invert the transformations) from the new input and update it to the output
        self.output_data_structure.UpdateData(self.PredictNeuralNetwork(data_structure_in = self.preprocessed_data_structure))
        output_value_list = np.squeeze(self.output_data_structure.ExportAsArray())

        # TODO: Right now, it only works with one variable that has the same shape as the output.
        if self.time >= self.time_buffer:
            output_value_index = 0
            for variable, source in self.dict_output.items():
            # Process related variables (e.g. TIME, STEP)
                if source == 'process':
                    self.output_model_part.ProcessInfo[variable] = output_value_list
                    output_value_index += 1
                # Node properties (e.g. position)
                elif source == 'node':
                    for node, node_id in zip(self.output_model_part.Nodes, range(self.output_model_part.NumberOfNodes())):
                        node.SetValue(variable, output_value_list)
                        output_value_index += 1
                # Node step values (e.g. variables like displacement)
                elif source == "solution_step":
                    for node, node_id in zip(self.output_model_part.Nodes, range(self.output_model_part.NumberOfNodes())):
                        node.SetSolutionStepValue(variable,0, output_value_list[node_id])
                        output_value_index += 1
                # Condition values
                elif source == "condition":
                    for condition, conditions_id in zip(self.output_model_part.GetConditions(), range(self.output_model_part.NumberOfConditions())):
                        condition.SetValue(variable, output_value_list[conditions_id])
                        output_value_index += 1

    def PredictNeuralNetwork(self, data_structure_in = None):
       
        data_out = NeuralNetworkData()
        
        self._PrintInfo("Predicting with the Neural Network...")

        if self.timesteps_as_features:
            data_structure_in.SetTimestepsAsFeatures()
        if self.feaures_as_timestpes:
            data_structure_in.SetFeaturesAsTimesteps()

        for process in self._GetListOfProcesses():
            try:
                output = process.Predict(self.neural_network_model, data_structure_in)
                if not output is None:
                    data_out.UpdateData(output)
            except AttributeError:
                pass
        if hasattr(data_structure_in, 'lookback_data'):
            data_structure_in.CheckLookbackAndUpdate(np.squeeze(data_out.ExportAsArray()))
        for process in self._GetListOfProcesses():
            try:
                prediction = process.TransformPredictions(self._GetListOfProcesses(), data_in = data_structure_in, data_out = data_out)
                if not prediction is None:
                    output_prediction = prediction
            except AttributeError:
                pass
        try:
            return output_prediction
        except UnboundLocalError:
            return data_out

    def FinalizeSolutionStep(self):

        if self.record:
            self.preprocessed_data_structure.data=self.record_data
        self.preprocessed_previous.UpdateData(self.preprocessed_data_structure.data)
        if not self.only_input:
            self.preprocessed_previous.UpdateLookbackAll(self.preprocessed_data_structure.lookback_data)
            self.preprocessed_previous.reorder_partitions = self.reorder_partitions
            self.reorder_partitions = 0 # Flag for only setting the reorder once
        self.converge_flag = True
    
    def AdvanceInTime(self, previous_time):
        self.time = previous_time + self.timestep
        if hasattr(self, 'model_geometry'):
            new_time = self.model_geometry.ProcessInfo[KratosMultiphysics.TIME] + self.timestep
            self.model_geometry.ProcessInfo.SetValue(KratosMultiphysics.TIME, new_time)
            new_step = self.model_geometry.ProcessInfo[KratosMultiphysics.STEP] + 1
            self.model_geometry.ProcessInfo.SetValue(KratosMultiphysics.STEP, new_step)
        return self.time

    def LoadProcessesList(self,list_of_processes):
        self._list_of_processes = list_of_processes

    def LoadNeuralNetworkModel(self, neural_network_model):
        self.neural_network_model = neural_network_model

    def LoadGeometry(self):
        self.model_geometry = self.model[self.model_geometry_name]

    def ImportModelPart(self):
        pass

    def _GetListOfProcesses(self):
        return self._list_of_processes
    
    def _PrintInfo(self, message):
        """This function prints info messages if the echo level is greater than 1.
        """
        if self.echo_level > 0:
            print(message)

    @staticmethod
    def _OrderSourcesFirst(values_list, variables_dictionary):
        """ Reorders the values list by sources instead of by variables (e.g. node by node)."""
        ordered_list = []
        k, m = divmod(len(values_list), len(variables_dictionary))
        # Split the lists
        split_list = list(values_list[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(len(variables_dictionary)))
        # Reorder the lists by source entry
        for index in range(len(split_list[0])):
            for variable_list in split_list:
                ordered_list.append(variable_list[index])

        return ordered_list

    @staticmethod
    def _OrderVariablesFirst(values_list, variables_dictionary):
        """ Reorders the values list by sources instead of by variables (e.g. node by node)."""
        ordered_list = []
        k, m = divmod(len(values_list), len(values_list)/len(variables_dictionary))
        # Split the lists
        split_list = list(values_list[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(len(values_list)/len(variables_dictionary)))
        # Reorder the lists by source entry
        for index in range(len(split_list[0])):
            for variable_list in split_list:
                ordered_list.append(variable_list[index])

        return ordered_list

    