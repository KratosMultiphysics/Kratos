from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

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
        if (type(model) != KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(project_parameters) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            self.is_printing_rank = (KratosMPI.mpi.rank == 0)
        else:
            self.is_printing_rank = True

        self._GetSolver().AddVariables() # this creates the solver and adds the variables

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.time < self.end_time:
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()


    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()

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
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self._GetSolver().Finalize()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()

        if self.is_printing_rank:
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
        # first we check if one of the output processes will print output in this step
        # this is done to save computation in case none of them will print
        is_output_step = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                is_output_step = True
                break

        if is_output_step: # at least one of the output processes will print output
            for process in self._GetListOfProcesses():
                process.ExecuteBeforeOutputStep()

            for output_process in self._GetListOfOutputProcesses():
                if output_process.IsOutputStep():
                    output_process.PrintOutput()

            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()


    def Check(self):
        """This function checks the AnalysisStage
        """
        pass

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

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")

    def _GetListOfProcesses(self):
        """This function returns the list of processes involved in this Analysis
        The list of processes is constructed in case it is not existing yet
        """
        if not hasattr(self, '_list_of_processes'):
            raise Exception("The list of processes was not yet created!")
        return self._list_of_processes

    def _GetListOfOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        The list of output processes is constructed in case it is not existing yet
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

        from process_factory import KratosProcessFactory
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
