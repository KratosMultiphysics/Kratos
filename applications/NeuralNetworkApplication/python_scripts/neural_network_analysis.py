from importlib import import_module
import KratosMultiphysics.NeuralNetworkApplication.input_dataclasses as InputDataclasses
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities
import KratosMultiphysics
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListDataWithLookback, ListNeuralNetworkData, NeuralNetworkData
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process_factory import NeuralNetworkProcessFactory

import tensorflow.keras as keras
import numpy as np

class NeuralNetworkAnalysis(AnalysisStage):
    """
    This class is the based on the AnalysisStage from Kratos but regarding the creation and training of the neural
    network model.
    """
    def __init__(self, project_parameters, model = None):      
        if not isinstance(project_parameters, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.project_parameters = project_parameters
        self.problem_type = self.project_parameters["problem_data"]["problem_type"].GetString()
        

        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        if self.problem_type == "predict_with_modelpart" or self.problem_type == "predict_without_modelpart":
            try:
                neural_network_file = self.project_parameters["problem_data"]["neural_network_file"].GetString()
                self._PrintInfo("Compiling neural network model from file...")
                self.neural_network_model = keras.models.load_model(neural_network_file)
                self.neural_network_model.compile()
            except AttributeError:
                self._PrintInfo("Model type predict must have a pre-trained model.")
        if self.problem_type == "predict_with_modelpart":
            try:
                self.solver_module = self.project_parameters["problem_data"]["solver_module"].GetString()
                self.solver_settings = self.project_parameters["problem_data"]["solver_settings"]
            except AttributeError:
                print("Warning: No solver especified.")
            try:
                self.model_geometry_file = self.project_parameters["problem_data"]["model_part_file"].GetString()
                self.model_geometry_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
                if model is None:
                    self.kratos_model = KratosMultiphysics.Model()
                    self.model_geometry = self.kratos_model.CreateModelPart(self.model_geometry_name)
                    KratosMultiphysics.ModelPartIO(self.model_geometry_file).ReadModelPart(self.model_geometry)
                    self.model_geometry.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = self.project_parameters["problem_data"]["solver_settings"]["solver_settings"]["domain_size"].GetInt()
                else:
                    self.kratos_model = model
                super().__init__(self.kratos_model, self.project_parameters)

            except AttributeError:
                raise Exception("The model part must be especified.")


            #TODO: This should be done through default parameters

            self.lookback = self.project_parameters["problem_data"] ["lookback"].GetInt()
            self.time_buffer = self.project_parameters["problem_data"]["time_buffer"].GetInt()
            self.timestep = self.project_parameters["problem_data"]["timestep"].GetDouble()
            self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()
            self.start_time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            try:
                self.record = self.project_parameters["problem_data"]["record"].GetBool()
            except RuntimeError:
                self.record = False
            try:
                self.only_input = self.project_parameters["problem_data"]["only_input"].GetBool()
            except RuntimeError:
                self.only_input = False
            try:
                self.timesteps_as_features = self.project_parameters["problem_data"]["timesteps_as_features"].GetBool()
            except RuntimeError:
                self.timesteps_as_features = False
            try:
                self.feaures_as_timestpes = self.project_parameters["problem_data"]["features_as_timesteps"].GetBool()
            except RuntimeError:
                self.feaures_as_timestpes = False
            try:
                self.soft_start_flag = self.project_parameters["problem_data"]["soft_start"].GetBool()
            except RuntimeError:
                self.soft_start_flag = True
            try:
                self.dimension_in = self.project_parameters["problem_data"]["dimension_input"].GetInt()
            except RuntimeError:
                self.dimension_in = 1
            try:
                self.reorder_partitions = self.project_parameters["problem_data"]["reorder_partitions"].GetInt()
            except RuntimeError:
                self.reorder_partitions = 1

            # TODO: This block is shared with data_generator_process, it could be separated and shared throug a function

            # getting the ModelPart from the Model
            output_model_part = self.project_parameters["problem_data"]["output_model_part"].GetString()
            if output_model_part == "":
                raise Exception('No "output_model_part" was specified!')
            self.output_model_part = self.kratos_model[output_model_part]

            # getting the input ModelPart from the Model
            input_model_part_name = self.project_parameters["problem_data"]["input_model_part"].GetString()
            if input_model_part_name == "":
                raise Exception('No "input_model_part" was specified!')
            self.input_model_part = self.kratos_model[input_model_part_name]

            # retrieving the input variables
            input_var_names = self.project_parameters["problem_data"]["input_variables"]
            variable_names = [ input_var_names[i].GetString() for i in range( input_var_names.size() ) ]
            input_sources_names = self.project_parameters["problem_data"]["input_sources"]
            self.input_sources = [ input_sources_names[i].GetString() for i in range( input_sources_names.size() ) ]
            self.input_variables = [ KratosMultiphysics.KratosGlobals.GetVariable( var ) for var in variable_names ]
            if len(self.input_variables) == 0:
                raise Exception('No variables specified for input!')
            if not (len(self.input_variables) == len(self.input_sources)):
                raise Exception('The number of input variables and sources are different.')
            self.dict_input = dict(zip(self.input_variables, self.input_sources))
            # getting input order
            try:
                self.input_order = self.project_parameters["problem_data"]["input_order"].GetString()
            except RuntimeError:
                self.input_order = "variables_first"

            # retrieving the output variables
            output_var_names = self.project_parameters["problem_data"]["output_variables"]
            variable_names = [ output_var_names[i].GetString() for i in range( output_var_names.size() ) ]
            output_sources_names = self.project_parameters["problem_data"]["output_sources"]
            self.output_sources = [ output_sources_names[i].GetString() for i in range( output_sources_names.size() ) ]
            self.output_variables = [ KratosMultiphysics.KratosGlobals.GetVariable( var ) for var in variable_names ]
            if len(self.output_variables) == 0:
                raise Exception('No variables specified for output!')
            if not (len(self.output_variables) == len(self.output_sources)):
                raise Exception('The number of output variables and sources are different.')
            self.dict_output = dict(zip(self.output_variables, self.output_sources))
            # getting output order
            try:
                self.output_order = self.project_parameters["problem_data"]["output_order"].GetString()
            except RuntimeError:
                self.output_order = "variables_first"

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        if self.problem_type == "predict_with_modelpart":
            super().Run()
        else:
            self.Initialize()
            self.RunTrainingLoop()
            self.Finalize()

    def Initialize(self, data_input = None, data_output = None):
        
        # Create list of processes
        self.__CreateListOfProcesses()
        
        if self.problem_type == "predict_with_modelpart":
            super().Initialize()
            self.input_data_structure = InputDataclasses.NeuralNetworkData()
            if self.lookback>0:
                self.preprocessed_data_structure = InputDataclasses.DataWithLookback(lookback_index=self.lookback)  
                self.preprocessed_previous = InputDataclasses.DataWithLookback(lookback_index=self.lookback)
            else:
                self.preprocessed_data_structure = InputDataclasses.NeuralNetworkData()
                self.preprocessed_previous = InputDataclasses.NeuralNetworkData()
            self.output_data_structure = InputDataclasses.NeuralNetworkData()
            self.time = 0.0

        else:
            if data_input != None:
                self.data_in = data_input
            else:
                self.data_in = []
            if data_input != None:
                self.data_out = data_output
            else:
                self.data_out = []

        # Preprocessing
        
            [self.data_in, self.data_out] = self.Preprocessing()
        
        ### TRAINING ###

        if self.problem_type == "train":

            self._PrintInfo("Starting generation of the neural network training setup...")

            # Initalize the network
            self._PrintInfo("Initializing the network...")
            inputs = self._GetListOfProcesses()[0].Initialize()
            outputs = inputs
            self._PrintInfo("Initialization done.")

            # Add layers to the network
            self._PrintInfo("Adding layers to the network...")
            for process in self._GetListOfProcesses()[1:]:
                outputs = process.Add(outputs)
            self._PrintInfo("Addition of the layers done.")

            # Generate the model and save it
            self._PrintInfo("Generating the model...")
            if not inputs == None:
                self.neural_network_model = keras.Model(inputs = inputs, outputs = outputs)
                self.neural_network_model.summary()
                for process in self._GetListOfProcesses():
                    process.Save(self.neural_network_model)
            else:
                self._PrintInfo("Warning: Model not defined.")
            self._PrintInfo("Generation of the model done.")

            # Compile the training parameters
            self._PrintInfo("Compiling the model(loss function and optimizers)...")
            self.loss_function = None
            self.optimizer = None
            for process in self._GetListOfProcesses():
                compilation_settings = process.Compile(self.loss_function, self.optimizer)
                if not compilation_settings == None:
                    self.loss_function = compilation_settings[0]
                    self.optimizer = compilation_settings[1]
            self._PrintInfo("Compilation of the model done.")

            # Compile the list of Callbacks
            self._PrintInfo("Compiling the callbacks...")
            self.callbacks_list = []
            for process in self._GetListOfProcesses():
                callback = process.Callback()
                if not callback == None:
                    self.callbacks_list.append(callback)
            self._PrintInfo("Compilation of the callbacks done.")

            # Compile the list of metrics
            self._PrintInfo("Compiling the metrics...")
            self.metrics_list = []
            for process in self._GetListOfProcesses():
                metric = process.CompileMetric()
                if not metric == None:
                    self.metrics_list.append(metric)
            self._PrintInfo("Compilation of the metrics done.")

            self._PrintInfo("Training model setup done.")

        ### TUNING ###

        elif self.problem_type == "tuning":

            self._PrintInfo("Starting generation of the neural network tuning setup...")
            def build_model(hp):
                # Initalize the network
                self._PrintInfo("Initializing the network...")
                inputs = self._GetListOfProcesses()[0].Initialize()
                outputs = inputs
                self._PrintInfo("Initialization done.")

                # Add layers to the network
                self._PrintInfo("Adding layers to the network...")
                for process in self._GetListOfProcesses()[1:]:
                    outputs = process.Add(outputs, hp = hp)
                self._PrintInfo("Addition of the layers done.")

                # Generate the model and save it
                self._PrintInfo("Generating the model...")
                model = keras.Model(inputs = inputs, outputs = outputs)
                self.loss_function = None
                self.optimizer = None
                for process in self._GetListOfProcesses():
                    compilation_settings = process.Compile(self.loss_function, self.optimizer, hp)
                    if not compilation_settings == None:
                        self.loss_function = compilation_settings[0]
                        self.optimizer = compilation_settings[1]
                self._PrintInfo("Generation of the model done.")

                # Compile the list of Callbacks
                self._PrintInfo("Compiling the callbacks...")
                self.callbacks_list = []
                for process in self._GetListOfProcesses():
                    callback = process.Callback()
                    if not callback == None:
                        self.callbacks_list.append(callback)
                self._PrintInfo("Compilation of the callbacks done.")
                
                # Compile the list of metrics
                self._PrintInfo("Compiling the metrics...")
                self.metrics_list = []
                for process in self._GetListOfProcesses():
                    metric = process.CompileMetric()
                    if not metric == None:
                        self.metrics_list.append(metric)
                                
                model.compile(optimizer = self.optimizer, loss = self.loss_function, metrics = self.metrics_list)
                self._PrintInfo("Compilation of the metrics done.")

                return model
            self.hypermodel = build_model
            self._PrintInfo("Tuning model setup done.")


        # Execute other processes that must run on the initialization
        self._PrintInfo("Executing processes that run on initialization...")
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        # Execute other processes that must run before the training
        if self.problem_type != "predict_with_modelpart":
            self._PrintInfo("Executing processes that run before the training...")
            for process in self._GetListOfProcesses():
                process.ExecuteBeforeTraining()

        # If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        self._PrintInfo("Executing processes that run on finalization...")
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self._PrintInfo("Transforming predictions if required...")
        for process in self._GetListOfProcesses():
            try:
                process.TransformPredictions(self._GetListOfProcesses())
            except AttributeError:
                pass
        
        self._PrintInfo("Plotting graphs if required...")
        for process in self._GetListOfProcesses():
            try:
                process.Plot()
            except AttributeError:
                pass

    def RunTrainingLoop(self):
        if self.problem_type == "train":
            self._PrintInfo("Training the model...")
            for process in self._GetListOfProcesses():
                process.ExecuteTraining(self.loss_function, self.optimizer, self.callbacks_list, self.metrics_list)
            self._PrintInfo("Training done.")
        elif self.problem_type == "tuning":
            self._PrintInfo("Tuning model...")
            for process in self._GetListOfProcesses():
                process.ExecuteTuning(self.hypermodel)
            self._PrintInfo("Tuning done.")

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.AdvanceInTime()
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def InitializeSolutionStep(self):
        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty

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
            
            # RUndo the reorder if indicated
        if self.output_order == 'sources_first':
            try:
                output_value_list = self._OrderVariablesFirst(output_value_list, self.dict_output.items())
            except IndexError:
                pass

        # Predict (and invert the transformations) from the new input and update it to the output
        self.output_data_structure.UpdateData(self.Predict(data_structure_in = self.preprocessed_data_structure))
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

    def FinalizeSolutionStep(self):

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

        if self.record:
            self.preprocessed_data_structure.data=self.record_data
        self.preprocessed_previous.UpdateData(self.preprocessed_data_structure.data)
        if not self.only_input:
            self.preprocessed_previous.UpdateLookbackAll(self.preprocessed_data_structure.lookback_data)
            self.preprocessed_previous.reorder_partitions = self.reorder_partitions
            self.reorder_partitions = 0 # Flag for only setting the reorder once
        self.converge_flag = True

    def Predict(self, data_structure_in = None):
       
        data_out = NeuralNetworkData()
        
        if self.problem_type in ["predict_with_modelpart","predict_without_modelpart", "predict"] :
            self._PrintInfo("Predicting with the Neural Network...")
            if data_structure_in is None:
                data_structure_in = self.data_in
            # if self.lookback >0:
            #     if not isinstance(data_structure_in, ListDataWithLookback):
            #         new_list = ListDataWithLookback()
            #         new_list.AddToList(data_structure_in)
            #         data_structure_in = new_list
            # else:
            #     if not isinstance(data_structure_in, ListNeuralNetworkData):
            #         new_list = ListNeuralNetworkData()
            #         new_list.AddToList(data_structure_in)
            #         data_structure_in = new_list
            
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
        else:
            raise Exception("The problem type for the neural network must be predict.")

    def Preprocessing(self, data_in = None, data_out = None):
        if not hasattr(self, '_list_of_processes'):
            self.__CreateListOfProcesses()
        if data_in is None:
            data_in = self.data_in
        if data_out is None:
            data_out = self.data_out
        self._PrintInfo("Starting data preprocessing...")
        if self.problem_type != "predict_with_modelpart":
            for process in self._GetListOfProcesses():
                preprocessing_settings = process.Preprocess(data_in, data_out)
                if not preprocessing_settings is None:
                    data_in = preprocessing_settings[0]
                    data_out = preprocessing_settings[1]
        self._PrintInfo("Data preprocessing done.")

        return [data_in, data_out]

    def AdvanceInTime(self):
        self.time += self.timestep
        if hasattr(self, 'model_geometry'):
            new_time = self.model_geometry.ProcessInfo[KratosMultiphysics.TIME] + self.timestep
            self.model_geometry.ProcessInfo.SetValue(KratosMultiphysics.TIME, new_time)
            new_step = self.model_geometry.ProcessInfo[KratosMultiphysics.STEP] + 1
            self.model_geometry.ProcessInfo.SetValue(KratosMultiphysics.STEP, new_step)


    def _CreateSolver(self):
        
        solver = import_module(self.solver_module).CreateSolver(self.kratos_model, self.solver_settings)
        return solver

    def _GetListOfProcesses(self):
        """This function returns the list of processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_processes'):
            raise Exception("The list of processes was not yet created!")
        return self._list_of_processes

    def _GetListOfOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_output_processes'):
            raise Exception("The list of output-processes was not yet created!")
        return self._list_of_output_processes

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of processes
        Format:
        "processes" : {
            initial_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ],
            boundary_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ]
        }
        The order of intialization can be specified by setting it in "initialization_order"
        if e.g. the "boundary_processes" should be constructed before the "initial_processes", then
        initialization_order should be a list containing ["boundary_processes", "initial_processes"]
        see the functions _GetOrderOfProcessesInitialization and _GetOrderOfOutputProcessesInitialization
        """
        list_of_processes = []

        factory = NeuralNetworkProcessFactory()

        if self.project_parameters.Has(parameter_name):
            processes_params = self.project_parameters[parameter_name]

            # first initialize the processes that depend on the order
            for processes_names in initialization_order:
                if processes_params.Has(processes_names):
                    list_of_processes += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in initialization_order:
                    if self.problem_type == "predict_with_modelpart":
                        list_of_processes += factory.ConstructListOfProcesses(value,self.kratos_model) # Does this work? or should it be processes[name]
                    else:
                        list_of_processes += factory.ConstructListOfProcesses(value) # Does this work? or should it be processes[name]

        return list_of_processes

    def _GetOrderOfProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the processes matters
        """
        return []

    def _GetOrderOfOutputProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the output-processes matters
        """
        return []

    def _GetSimulationName(self):
        """Returns the name of the Simulation
        """
        return "Analysis"

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes

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
