from importlib import import_module
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities
import KratosMultiphysics
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import  NeuralNetworkData
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process_factory import NeuralNetworkProcessFactory
from KratosMultiphysics.NeuralNetworkApplication import python_solvers_wrapper_neural_network as neural_network_solvers
from KratosMultiphysics.NeuralNetworkApplication.neural_network_model import MachineLearningModel
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
        try:
            self.ml_model = MachineLearningModel(self.project_parameters["problem_data"])
        except RuntimeError:
            self.ml_model = MachineLearningModel('keras')

        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        if self.project_parameters["problem_data"].Has("mesh_moving"):
            self.mesh_moving = self.project_parameters["problem_data"]["mesh_moving"].GetBool()
        else:
            self.mesh_moving = False

        if self.problem_type == "predict_with_modelpart" or self.problem_type == "predict_without_modelpart":
            try:
                neural_network_file = self.project_parameters["problem_data"]["neural_network_file"].GetString()
                self._PrintInfo("Compiling neural network model from file...")
                self.neural_network_model = self.ml_model.load_model(neural_network_file)
            except AttributeError:
                self._PrintInfo("Model type predict must have a pre-trained model.")
        if self.problem_type == "predict_with_modelpart":
            try:
                self.solver_module = self.project_parameters["problem_data"]["solver_module"].GetString()
                self.solver_settings = self.project_parameters["problem_data"]["solver_settings"]
            except AttributeError:
                print("Warning: No solver especified.")
            # try:
            self.model_geometry_file = self.project_parameters["problem_data"]["model_part_file"].GetString()
            self.model_geometry_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
            if model is None:
                self.kratos_model = KratosMultiphysics.Model()
                self.model_geometry = self.kratos_model.CreateModelPart(self.model_geometry_name)
                self.solver = self._CreateNeuralNetworkSolver()
                super().__init__(self.kratos_model, self.project_parameters)
                if self.mesh_moving:
                    self._GetNeuralNetworkSolver().AddVariables()
                # KratosMultiphysics.ModelPartIO(self.model_geometry_file).ReadModelPart(self.model_geometry)
                try:
                    self.model_geometry.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = self.project_parameters["problem_data"]["solver_settings"]["solver_settings"]["domain_size"].GetInt()
                except RuntimeError:
                    self.model_geometry.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = self.project_parameters["problem_data"]["domain_size"].GetInt()
            else:
                self.kratos_model = model
                self.model_geometry = self.kratos_model.CreateModelPart(self.model_geometry_name)
                self.solver = self._CreateNeuralNetworkSolver()
                super().__init__(self.kratos_model, self.project_parameters)
                if self.mesh_moving:
                    self._GetNeuralNetworkSolver().AddVariables()
                # KratosMultiphysics.ModelPartIO(self.model_geometry_file).ReadModelPart(self.kratos_model[self.model_geometry_name])
                try:
                    self.model_geometry.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = self.project_parameters["problem_data"]["solver_settings"]["solver_settings"]["domain_size"].GetInt()
                except RuntimeError:
                    self.model_geometry.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = self.project_parameters["problem_data"]["domain_size"].GetInt()
            

    def _CreateSolver(self):
        """This function creates the PHYSICS solver."""
        return import_module(self.solver_module).CreateSolver(self.kratos_model, self.solver_settings)

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
        
        if self.problem_type == "predict_with_modelpart":
            print(self.kratos_model)
            self._CreateModelers()
            self._ModelersSetupGeometryModel()
            self._ModelersPrepareGeometryModel()
            self._ModelersSetupModelPart()

            self._GetSolver().ImportModelPart()
            self._GetSolver().PrepareModelPart()
            if self.mesh_moving:
                self._GetNeuralNetworkSolver().AddDofs()
            self._GetSolver().AddDofs()

            self.ModifyInitialProperties()
            self.ModifyInitialGeometry()

            ##here we initialize user-provided processes
            self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
            for process in self._GetListOfProcesses():
                process.ExecuteInitialize()

            self._GetSolver().Initialize()
            self._GetNeuralNetworkSolver().Initialize()
            self.Check()

            self.ModifyAfterSolverInitialize()

            for process in self._GetListOfProcesses():
                process.ExecuteBeforeSolutionLoop()

            ## Stepping and time settings
            self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

            if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            else:
                self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

            ## If the echo level is high enough, print the complete list of settings used to run the simualtion
            if self.echo_level > 1:
                with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                    parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")
            
        else:
            if data_input != None:
                self.data_in = data_input
            else:
                self.data_in = []
            if data_output != None:
                self.data_out = data_output
            else:
                self.data_out = []

        # Create list of processes

        # Preprocessing
        
            [self.data_in, self.data_out] = self.Preprocessing()
        
            self.__CreateListOfProcesses()

        if self.problem_type == "predict_with_modelpart":
            self.solver.LoadProcessesList(self._GetListOfProcesses())
            self.solver.LoadNeuralNetworkModel(self.neural_network_model)

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
                self.neural_network_model = self.ml_model.create_model(inputs = inputs, outputs = outputs)
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
                model = self.ml_model.create_model(inputs = inputs, outputs = outputs)
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
        if not self.problem_type == "predict_with_modelpart":
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
            self.time = self._GetNeuralNetworkSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetNeuralNetworkSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def InitializeSolutionStep(self):
        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty

        self._GetNeuralNetworkSolver().InitializeSolutionStep()
       

    def FinalizeSolutionStep(self):
        
        self._GetNeuralNetworkSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

    def Predict(self, data_structure_in = None):
       
       # TODO: This is no longer neccessary, it is only for legacy wrappers. It should be deleted.

        data_out = NeuralNetworkData()
        
        if self.problem_type in ["predict_with_modelpart","predict_without_modelpart", "predict"] :
            self._PrintInfo("Predicting with the Neural Network...")
            if data_structure_in is None:
                data_structure_in = self.data_in
            

            for process in self._GetListOfProcesses():
                try:
                    output = process.Predict(self.ml_model.model, data_structure_in)
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

    def _CreateNeuralNetworkSolver(self):
        return neural_network_solvers.CreateSolver(self.kratos_model, self.project_parameters["problem_data"])


    def _GetNeuralNetworkSolver(self):
        """This function gets the NEURAL NETWORK solver."""
        return self.solver

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
        return "Neural Network Analysis"

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

    
