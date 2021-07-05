import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process_factory import NeuralNetworkProcessFactory

import tensorflow.keras as keras

class NeuralNetworkAnalysis(AnalysisStage):
    """
    This class is the based on the AnalysisStage from Kratos but regarding the creation and training of the neural
    network model.
    """
    def __init__(self, project_parameters):      
        if not isinstance(project_parameters, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.project_parameters = project_parameters
        self.problem_type = self.project_parameters["problem_data"]["problem_type"].GetString()
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunTrainingLoop()
        self.Finalize()

    def Initialize(self):
        # Create list of processes
        self.__CreateListOfProcesses()

        # Preprocessing
        self.data_in = None
        self.data_out = None
        for process in self._GetListOfProcesses():
            preprocessing_settings = process.Preprocess(self.data_in, self.data_out)
            if not preprocessing_settings is None:
                self.data_in = preprocessing_settings[0]
                self.data_out = preprocessing_settings[1]

        if self.problem_type == "train":
            # Initalize the network
            inputs = self._GetListOfProcesses()[0].Initialize()
            outputs = inputs

            # Add layers to the network
            for process in self._GetListOfProcesses()[1:]:
                outputs = process.Add(outputs)

            # Generate the model and save it
            if not inputs == None:
                self.model = keras.Model(inputs = inputs, outputs = outputs)
                self.model.summary()
                for process in self._GetListOfProcesses():
                    process.Save(self.model)
            else:
                print("Warning: Model not defined.")

            # Compile the training parameters
            self.loss_function = None
            self.optimizer = None
            for process in self._GetListOfProcesses():
                compilation_settings = process.Compile(self.loss_function, self.optimizer)
                if not compilation_settings == None:
                    self.loss_function = compilation_settings[0]
                    self.optimizer = compilation_settings[1]

            # Compile the list of Callbacks
            self.callbacks_list = []
            for process in self._GetListOfProcesses():
                callback = process.Callback()
                if not callback == None:
                    self.callbacks_list.append(callback)
            
            # Compile the list of metrics
            self.metrics_list = []
            for process in self._GetListOfProcesses():
                metric = process.CompileMetric()
                if not metric == None:
                    self.metrics_list.append(metric)

        elif self.problem_type == "tuning":
            def build_model(hp):
                # Initalize the network
                inputs = self._GetListOfProcesses()[0].Initialize()
                outputs = inputs

                # Add layers to the network
                for process in self._GetListOfProcesses()[1:]:
                    outputs = process.Add(outputs, hp = hp)

                # Generate the model and save it
                model = keras.Model(inputs = inputs, outputs = outputs)
                self.loss_function = None
                self.optimizer = None
                for process in self._GetListOfProcesses():
                    compilation_settings = process.Compile(self.loss_function, self.optimizer, hp)
                    if not compilation_settings == None:
                        self.loss_function = compilation_settings[0]
                        self.optimizer = compilation_settings[1]

                # Compile the list of Callbacks
                self.callbacks_list = []
                for process in self._GetListOfProcesses():
                    callback = process.Callback()
                    if not callback == None:
                        self.callbacks_list.append(callback)
                
                # Compile the list of metrics
                self.metrics_list = []
                for process in self._GetListOfProcesses():
                    metric = process.CompileMetric()
                    if not metric == None:
                        self.metrics_list.append(metric)
                
                model.compile(optimizer = self.optimizer, loss = self.loss_function, metrics = self.metrics_list)

                return model
            self.hypermodel = build_model


        # Execute other processes that must run on the initialization
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        # Execute other processes that must run before the training
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
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()
        for process in self._GetListOfProcesses():
            process.TransformPredictions(self._GetListOfProcesses())
        for process in self._GetListOfProcesses():
            process.Plot()

    def RunTrainingLoop(self):
        if self.problem_type == "train":
            for process in self._GetListOfProcesses():
                process.ExecuteTraining(self.loss_function, self.optimizer, self.callbacks_list, self.metrics_list)
        elif self.problem_type == "tuning":
            for process in self._GetListOfProcesses():
                process.ExecuteTuning(self.hypermodel)

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
