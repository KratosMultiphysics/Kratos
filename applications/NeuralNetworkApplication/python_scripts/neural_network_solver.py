from importlib import import_module
from KratosMultiphysics import python_solver
import KratosMultiphysics.NeuralNetworkApplication.input_dataclasses as InputDataclasses
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities
import KratosMultiphysics
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListDataWithLookback, ListNeuralNetworkData, NeuralNetworkData
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process_factory import NeuralNetworkProcessFactory
from importlib import import_module

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
        self.lookback = project_parameters ["lookback"].GetInt()
        self.time_buffer = project_parameters["time_buffer"].GetDouble()
        self.timestep = project_parameters["timestep"].GetDouble()
        self.end_time = project_parameters["end_time"].GetDouble()
        self.start_time = project_parameters["start_time"].GetDouble()

        if project_parameters.Has("import_applications"):
            applications_list = project_parameters["import_applications"].GetStringArray()
            for application in applications_list:
                import_module("KratosMultiphysics." + application)
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

        # TODO: This block is shared with data_generator_process, it could be separated and shared throug a function

        # getting the ModelPart from the Model
        self.output_model_parts = []
        if not project_parameters.Has("output_model_part"):
            raise Exception("Output model part not present")
        if project_parameters["output_model_part"].IsArray():
            self.output_is_vector = True
            self.output_model_parts_names = project_parameters["output_model_part"].GetStringArray()

        else:
            self.output_is_vector = False
            self.output_model_parts_names = [project_parameters["output_model_part"].GetString()]
        if not self.output_model_parts_names:
            raise Exception('No "output_model_part" was specified!')

        # getting the input ModelPart from the Model
        self.input_model_parts = []
        if not project_parameters.Has("input_model_part"):
            raise Exception("Input model part not present")
        if project_parameters["input_model_part"].IsArray():
            self.input_is_vector = True
            self.input_model_parts_names = project_parameters["input_model_part"].GetStringArray()
        else:
            self.input_is_vector = False
            self.input_model_parts_names = [project_parameters["input_model_part"].GetString()]
        if not self.input_model_parts_names:
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
        for name in self.output_model_parts_names:
            self.output_model_parts.append(self.model[name])
        if not (len(self.output_variables) == len(self.output_model_parts)):
            if self.output_is_vector:
                raise Exception('The output_model_parts are given as vector -- The number of output variables and model parts must be the same.')
            else:
                while not (len(self.output_variables) == len(self.output_model_parts)):
                    self.output_model_parts.append(self.output_model_parts[-1])
        self.dict_output = list(zip(self.output_model_parts, self.output_variables, self.output_sources))
        
        for name in self.input_model_parts_names:
            self.input_model_parts.append(self.model[name])
        if not (len(self.input_variables) == len(self.input_model_parts)):
            if self.input_is_vector:
                raise Exception('The input_model_parts are given as vector -- The number of input variables and model parts must be the same.')
            else:
                while not (len(self.input_variables) == len(self.input_model_parts)):
                    self.input_model_parts.append(self.input_model_parts[-1])
        self.dict_input = list(zip(self.input_model_parts, self.input_variables, self.input_sources))



    def InitializeSolutionStep(self):

        pass


    def SolveSolutionStep(self):

        print("Receiving data into the neural network model")
        current_time = self.input_model_parts[0].ProcessInfo[KratosMultiphysics.TIME]
        input_value_list=[]
        # TODO: Modify to use vector variables instead of the components (unzip Kratos.Array)

        previous_timestep_model_part_index = 0
        for model_part, variable, source in self.dict_input:
            model_input_value_list = []
            # Process related variables (e.g. TIME, STEP)
            if source == 'process':
                model_input_value_list.append(model_part.ProcessInfo[variable])
            # Node properties (e.g. position)
            elif source == 'node':
                for node in model_part.Nodes:
                    model_input_value_list.append(getattr(node,variable.Name()))
            # Node step values (e.g. variables like displacement)
            elif source == "solution_step":
                for node in model_part.Nodes:
                    model_input_value_list.append(node.GetSolutionStepValue(variable,0))
            elif source == "previous_solution_step":
                # TODO: modify to enable the use of the surrogate model only in a set interval
                # and the physics-based one in the rest. Other parts of the code are also affected.
                if current_time > self.timestep:
                    model_input_value_list = self.input_from_previous[previous_timestep_model_part_index]
                    previous_timestep_model_part_index += 1
                else:
                    for node in model_part.Nodes:
                        model_input_value_list.append(node.GetSolutionStepValue(variable,0))
            # Condition values
            elif source == "condition":
                for condition in model_part.GetConditions():
                    model_input_value_list.append(condition.GetValue(variable))
            input_value_list.extend(model_input_value_list)
        # Reorder if indicated
        if self.input_order == 'sources_first':
            try:
                model_input_value_list = self._OrderSourcesFirst(model_input_value_list, self.dict_input)
            except IndexError:
                pass
        self.input_from_modelpart = np.array(input_value_list)

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
            output_value_list = []
            for model_part, variable, source in self.dict_output:
                output_model_value_list = []
                # Process related variables (e.g. TIME, STEP)
                if source == 'process':
                    output_model_value_list.append(model_part.ProcessInfo[variable])
                # Node properties (e.g. position)
                elif source == 'node':
                    for node in model_part.Nodes:
                        output_model_value_list.append(getattr(node,variable.Name()))
                # Node step values (e.g. variables like displacement)
                elif source == "solution_step":
                    for node in model_part.Nodes:
                        output_model_value_list.append(node.GetSolutionStepValue(variable,0))
                elif source == "previous_solution_step":
                    for node in model_part.Nodes:
                        output_model_value_list.append(node.GetSolutionStepValue(variable,0))
                # Condition values
                elif source == "condition":
                    for condition in model_part.GetConditions():
                        output_model_value_list.append(condition.GetValue(variable))
                output_value_list.extend(output_model_value_list)
            # Reorder if indicated
            if self.output_order == 'sources_first':
                try:
                    output_value_list = self._OrderSourcesFirst(output_value_list, self.dict_output)
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

            self.preprocessed_data_structure.UpdateData(np.transpose(self.input_data_structure.ExportAsArray()[:,np.newaxis]))
            # self.preprocessed_data_structure.UpdateData(self.input_data_structure.ExportAsArray())
        # Set the reordering of the preprocessed data structure
        if self.reorder_partitions > 1:
            self.preprocessed_data_structure.reorder_partitions = self.reorder_partitions
            

        # Predict (and invert the transformations) from the new input and update it to the output
        self.output_data_structure.UpdateData(self.PredictNeuralNetwork(data_structure_in = self.preprocessed_data_structure))
        output_value_list = np.squeeze(self.output_data_structure.ExportAsArray())

            # Undo the reorder if indicated
        if self.output_order == 'sources_first':
            try:
                output_value_list = self._OrderVariablesFirst(output_value_list, self.dict_output)
            except IndexError:
                pass
        # TODO: Right now, it only works with one variable that has the same shape as the output.
        if self.time >= self.time_buffer: # TODO: Check for consistency with timesteps nad time
            output_value_index = 0
            for model_part, variable, source in self.dict_output:
            # Process related variables (e.g. TIME, STEP)
                if source == 'process':
                    model_part.ProcessInfo[variable] = output_value_list
                    output_value_index += 1
                # Node properties (e.g. position)
                elif source == 'node':
                    for node, node_id in zip(model_part.Nodes, range(model_part.NumberOfNodes())):
                        node.SetValue(variable, output_value_list[node_id + output_value_index])
                    output_value_index += node_id
                # Node step values (e.g. variables like displacement)
                elif source == "solution_step":
                    for node, node_id in zip(model_part.Nodes, range(model_part.NumberOfNodes())):
                        node.SetSolutionStepValue(variable,0, output_value_list[node_id + output_value_index])
                    output_value_index += node_id
                # Condition values (should be checked)
                elif source == "condition":
                    for condition, conditions_id in zip(model_part.GetConditions(), range(model_part.NumberOfConditions())):
                        condition.SetValue(variable, output_value_list[conditions_id + output_value_index])
                    output_value_index += conditions_id

    def PredictNeuralNetwork(self, data_structure_in = None):
       
        data_out = NeuralNetworkData()
        
        self._PrintInfo("Predicting with the Neural Network...")

        if self.timesteps_as_features:
            data_structure_in.SetTimestepsAsFeatures()
        if self.feaures_as_timestpes:
            data_structure_in.SetFeaturesAsTimesteps()

        if self.external_model:
        # NOTE: This is needed due to unexpected bug with the predict process.
        # In normal conditions, it is not necessary, but it must be used if there are 
        # processes that parse from a function in string form.
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

    def Preprocessing(self, data_in = None, data_out = None):
        if not hasattr(self, '_list_of_processes'):
            self.__CreateListOfProcesses()
        self._PrintInfo("Starting data preprocessing...")
        for process in self._GetListOfProcesses():
            try:
                preprocessing_settings = process.Preprocess(data_in, data_out)
            except AttributeError:
                preprocessing_settings = None

            if not preprocessing_settings is None:
                data_in = preprocessing_settings[0]
                data_out = preprocessing_settings[1]
        self._PrintInfo("Data preprocessing done.")

        return [data_in, data_out]

    def FinalizeSolutionStep(self):

        if self.record:
            self.preprocessed_data_structure.data=self.record_data
        self.preprocessed_previous.UpdateData(self.preprocessed_data_structure.data)
        if not self.only_input:
            self.preprocessed_previous.UpdateLookbackAll(self.preprocessed_data_structure.lookback_data)
            self.preprocessed_previous.reorder_partitions = self.reorder_partitions
            self.reorder_partitions = 0 # Flag for only setting the reorder once

        self.input_from_previous = []
        for model_part, variable, source in self.dict_input:
            model_input_value_list = []
            if source == "previous_solution_step":
                for node in model_part.Nodes:
                    model_input_value_list.append(node.GetSolutionStepValue(variable,0))
                self.input_from_previous.append(model_input_value_list)
                
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
        if self.echo_level > 0:
            self.neural_network_model.summary()

    def LoadGeometry(self):
        self.model_geometry = self.model[self.model_geometry_name]

    def ImportModelPart(self):
        self._ImportModelPart(self.model[self.model_geometry_name], self.model_import_settings)

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
        k, m = int(k), int(m)
        # Split the lists
        split_list = list(values_list[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(int(len(values_list)/len(variables_dictionary))))
        # Reorder the lists by source entry
        for index in range(len(split_list[0])):
            for variable_list in split_list:
                ordered_list.append(variable_list[index])

        return ordered_list

    