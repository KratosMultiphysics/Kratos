# Importing Python modules
from importlib import import_module

# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

class AnalysisStage(object):
    """The base class for the AnalysisStage-classes in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, project_parameters):
        """The constructor of the AnalysisStage-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        project_parameters -- The ProjectParameters used
        """
        if not isinstance(model, KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if not isinstance(project_parameters, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()
        is_distributed_run = KratosMultiphysics.IsDistributedRun()

        if self.parallel_type == "OpenMP" and is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"OpenMP" is specified as "parallel_type", but Kratos is running distributed!')
        if self.parallel_type == "MPI" and not is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"MPI" is specified as "parallel_type", but Kratos is not running distributed!')

        self._GetSolver().AddVariables() # this creates the solver and adds the variables
        self.__AddProcessesVariables() # this adds the variables required by the processes

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.time < self.end_time

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()

        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()
        self.__AddProcessesDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
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
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time

        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self._GetSolver().Finalize()

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")

    def GetFinalData(self):
        """Returns the final data dictionary.

        The main purpose of this function is to retrieve any data (in a key-value format) from outside the stage.
        Note that even though it can be called at any point, it is intended to be called at the end of the stage run.
        """

        return {}

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()


    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        execute_was_called = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                if not execute_was_called:
                    for process in self._GetListOfProcesses():
                        process.ExecuteBeforeOutputStep()
                    execute_was_called = True

                output_process.PrintOutput()

        if execute_was_called:
            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

    def Check(self):
        """This function checks the AnalysisStage

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Checking solver
        self._GetSolver().Check()

        # Checking processes
        for process in self._GetListOfProcesses():
            process.Check()

    def Clear(self):
        """This function clears the AnalysisStage

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Clearing solver
        self._GetSolver().Clear()

        # Clearing processes
        for process in self._GetListOfProcesses():
            process.Clear()

    def ModifyInitialProperties(self):
        """this is the place to eventually modify material properties in the stage """
        pass

    def ModifyInitialGeometry(self):
        """this is the place to eventually modify geometry (for example moving nodes) in the stage """
        pass

    def ModifyAfterSolverInitialize(self):
        """this is the place to eventually do any modification that requires the solver to be initialized """
        pass

    def ApplyBoundaryConditions(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        #other operations as needed

    def ChangeMaterialProperties(self):
        """this function is where the user could change material parameters as a part of the solution step """
        pass

    def Save(self, serializer: KratosMultiphysics.StreamSerializer) -> None:
        """Serializes current analysis stage instance

        This method is intended to make the class pure Python (pickable). This means serialize all the Kratos objects,
        that is to say all the objects coming from Pybind, with the provided serializer. After the serialization, it is
        required to assign None value to all the objects in order to make the class pickable.
        """
        pass

    def Load(self, serializer: KratosMultiphysics.StreamSerializer) -> None:
        """Loads current analysis stage instance

        From the given serializer, this method restores current class from a pure Python status (pickable) to the one in the serializer.
        """
        pass

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")

    def _AdvanceTime(self):
        """ Computes the following time
            The default method simply calls the solver
        """
        return self._GetSolver().AdvanceInTime(self.time)

    ### Modelers
    def _ModelersSetupGeometryModel(self):
        # Import or generate geometry models from external input.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model started.")
            modeler.SetupGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model finished.")

    def _ModelersPrepareGeometryModel(self):
        # Prepare or update the geometry model_part.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model started.")
            modeler.PrepareGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model finished.")

    def _ModelersSetupModelPart(self):
        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

    ### Modelers
    def _GetListOfModelers(self):
        """ This function returns the list of modelers
        """
        if not hasattr(self, '_list_of_modelers'):
            raise Exception("The list of modelers was not yet created!")
        return self._list_of_modelers

    def _CreateModelers(self):
        """ List of modelers in following format:
        "modelers" : [{
            "name" : "geometry_import",
            "parameters" : {
                "echo_level" : 0,
                // settings for this modeler
            }
        },{ ... }]
        """
        self._list_of_modelers = []

        if self.project_parameters.Has("modelers"):
            modelers_list = self.project_parameters["modelers"]
            #TODO: Remove this after the deprecation period
            if self.__BackwardCompatibleModelersCreation(modelers_list):
                from KratosMultiphysics.modeler_factory import KratosModelerFactory
                factory = KratosModelerFactory()
                self._list_of_modelers = factory.ConstructListOfModelers(self.model, modelers_list)
            else:
                factory = KratosModelParametersFactory(self.model)
                self._list_of_modelers = factory.ConstructListOfItems(modelers_list)

    @classmethod
    def __BackwardCompatibleModelersCreation(self, modelers_list):
        return any([modeler.Has("modeler_name") for modeler in modelers_list])

    ### Processes
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
        The order of initialization can be specified by setting it in "initialization_order"
        if e.g. the "boundary_processes" should be constructed before the "initial_processes", then
        initialization_order should be a list containing ["boundary_processes", "initial_processes"]
        see the functions _GetOrderOfProcessesInitialization and _GetOrderOfOutputProcessesInitialization
        """
        list_of_processes = []

        factory = KratosProcessFactory(self.model)

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

    def _CheckDeprecatedOutputProcesses(self, list_of_processes):
        deprecated_output_processes = []
        for process in list_of_processes:
            if issubclass(type(process), KratosMultiphysics.OutputProcess):
                deprecated_output_processes.append(process)
                msg  = "{} is an OutputProcess. However, it has been constructed as a regular process.\n"
                msg += "Please, define it as an 'output_processes' in the ProjectParameters."
                IssueDeprecationWarning("AnalysisStage", msg.format(process.__class__.__name__))
        return deprecated_output_processes

    def _GetSimulationName(self):
        """Returns the name of the Simulation
        """
        return "Analysis"

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

    def __AddProcessesVariables(self):
        variables_list = []
        if self.project_parameters.Has("processes"): #TODO: Here we are assuming that output_processes cannot add variables. Decide this in KTC.
            # Get the standard (not output) processes
            processes_settings = self.project_parameters["processes"]

            # Loop each list of processes (e.g., initial_conditions, boundary_conditions, etc.).
            for _, processes_list in processes_settings.items():
                for process in processes_list:
                    # Get the process input settings
                    # Note that we support both 'parameters' and 'Parameters' to keep processes backward compatibility
                    if process.Has("parameters"):
                        parameters = process["parameters"]
                    elif process.Has("Parameters"):
                        IssueDeprecationWarning("AnalysisStage", f"Found 'Parameters' in {process} definition. Use 'parameters' instead.")
                        parameters = process["Parameters"]
                    else:
                        raise NameError(f"Provided item '{process}' is incomplete.'parameters' field must be defined for each item.")

                    # Get the process class
                    if process.Has("name"):
                        # Check if current process is registered
                        registry_entry = process["name"].GetString()
                        if KratosMultiphysics.Registry.HasItem(registry_entry):
                            # Get class from stored Python module
                            if KratosMultiphysics.Registry.HasItem(f"{registry_entry}.ModuleName"):
                                class_name = registry_entry.split(".")[-1]
                                module_name = KratosMultiphysics.Registry[f"{registry_entry}.ModuleName"]
                                module = import_module(module_name)

                                # Check if the standard class name is available
                                # Note that in here we assume that the registry last key matches the class name
                                if hasattr(module, class_name):
                                    class_attribute = getattr(module, class_name)
                                else:
                                    # If the registry key does not contain the class name we check for the 'ClassName' entry
                                    if KratosMultiphysics.Registry.HasItem(f"{registry_entry}.ClassName"):
                                        class_name = KratosMultiphysics.Registry[f"{registry_entry}.ClassName"]
                                        KratosMultiphysics.Logger.PrintWarning(f"Trying to get item from non-standard 'ClassName' value {class_name}.")
                                        class_attribute = getattr(module, class_name)
                                    else:
                                        err_msg = f"The '{class_name}' class name cannot be found within the '{module_name}' module."
                                        raise Exception(err_msg)
                            else:
                                raise Exception(f"Registry process '{registry_entry}' cannot be obtained.")

                            # Call the static specifications method of the obtained class
                            if hasattr(class_attribute, "GetSpecifications"):
                                specifications = class_attribute.GetSpecifications(parameters)
                                variables_list.extend(specifications["required_solution_step_data_variables"].GetStringArray())
                        else:
                            KratosMultiphysics.Logger.PrintWarning(f"Asking to retrieve the non-registered 'name' '{registry_entry}'.")
                    else:
                        # Alternative (old) retrieval from Python module
                        if not process.Has("python_module"):
                            raise NameError(f'"python_module" must be defined in your parameters. Check all your processes')

                        python_module_name = process["python_module"].GetString() # python-script that contains the process
                        if process.Has("kratos_module"): # for Kratos-processes
                            # Get the class name from the Python module
                            kratos_module_name = process["kratos_module"].GetString()
                            if not kratos_module_name.startswith("KratosMultiphysics"):
                                kratos_module_name = "KratosMultiphysics." + kratos_module_name
                            module_name = kratos_module_name + "." + python_module_name
                        else: # for user-defined processes
                            module_name = import_module(python_module_name)

                        # Obtain the class attribute
                        # Note that in here we assume that the Python module name is the class name in snake case
                        class_name = ''.join(word.title() for word in python_module_name.split('_'))
                        module = import_module(module_name)
                        class_attribute = getattr(module, class_name)

                        # Call the static specifications method of the obtained class
                        # Retrieve the required variables list from the specifications and append them to the list
                        if hasattr(class_attribute, "GetSpecifications"):
                            specifications = class_attribute.GetSpecifications(parameters)
                            variables_list.extend(specifications["required_solution_step_data_variables"].GetStringArray())

        # Make the list unique
        variables_list = list(set(variables_list))

        # Add the variables to the solver computing model part
        for variable_name in variables_list:
            variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
            self._GetSolver().GetComputingModelPart().AddNodalSolutionStepVariable(variable)

    def __AddProcessesDofs(self):
        pass

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning

        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("AnalysisStage", warn_msg)
