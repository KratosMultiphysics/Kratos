# co simulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.base_co_simulation_classes.co_simulation_base_solver import CoSimulationBaseSolver
# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_data_structure as data_str
cs_data_structure = data_str.__DATA_STRUCTURE__
import collections

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
#  This class is intended to server as the base class for all the coupled solvers.
class CoSimulationBaseCoupledSolver(CoSimulationBaseSolver):
    ## The constructor
    #

    #  @param custom_settings     parameters for configuring the CoSimulationBaseCoupledSolver
    def __init__(self, custom_settings):
        ##settings string in json format
        # default values for all available settings
        # for mandatory settings, the type is defined
        self.full_settings = custom_settings
        self.settings= custom_settings['coupled_solver_settings']
        #super(CoSimulationBaseCoupledSolver,self).__init__(custom_settings)
        default_setting = cs_data_structure.Parameters("""
        {
            "name" : "",
            "solver_type" : "gauss_seidel_strong_coupling",
            "echo_level" : 0,
            "num_coupling_iterations" : 10,
            "start_coupling_time" : 0.0,
            "convergence_accelerators" : [],
            "participants" : [],
            "convergence_criteria_settings" : {}
        }
        """)
        #Comment predictors are missing (I think those can be optional, but also used for weak-coupling)
        #Comment settings for predictor, conv_crit and accelerator should be consitently defined!
        self.settings.ValidateAndAssignDefaults(default_setting)
        self.number_of_participants = self.settings['participants'].size()
        self.echo_level = self.settings["echo_level"].GetInt()
        self.num_coupling_iterations = self.settings["num_coupling_iterations"].GetInt() #Comment should not be here

        # Get the participating solvers a map with their names and objects
        self.participating_solvers = self._CreateSolvers(self.full_settings['solvers'])
        self.solver_settings = self._GetSolverCoSimulationDetails(self.settings["participants"])

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if self.settings.Has("start_coupling_time"):
           self.start_coupling_time = self.settings["start_coupling_time"].GetDouble()
        if self.start_coupling_time > 0.0:
            self.coupling_started = False
        else:
            self.coupling_started = True
    ## Initialize : Initialize function. Called only once
    #               all member variables are initialized here.

    def Initialize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Initialize()

    ## Finalize : Finalize function. Called only once
    #               all member variables are finalized here.

    def Finalize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Finalize()

    ## Predict : Predict the solution of the next solution step.
    #

    def Predict(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Predict()

    ## InitializeSolutionStep : Called once in the beginning of the solution step
    #

    def InitializeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.InitializeSolutionStep()

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #

    def FinalizeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.FinalizeSolutionStep()

    ## OutputSolutionStep : Called once at the end of the solution step.
    #                       The output of the solvers and / or co-simulation output
    #                       can be performed in this function.
    #

    def OutputSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.OutputSolutionStep()

    ## AdvanceInTime :  This function will advance the current solver to the given current_time
    #
    #  @param self                      The object pointer.
    #  @current_time                    The time to which the current solver should advance
    def AdvanceInTime(self, current_time):
        new_time = 0.0
        for solver_name, solver in self.participating_solvers.items():
            new_time = max(solver.AdvanceInTime(current_time), new_time)

        if self.start_coupling_time > new_time:
            self.coupling_started = False
        else:
            self.coupling_started = True

        return new_time

    ## SolveSolutionStep : This function implements the coupling workflow
    #                      specific to the coupled solver by using
    #                      functions of its participating solvers in a specific
    #                      order.
    #
    def SolveSolutionStep(self):
        err_msg  = 'Calling "SolveSolutionStep" of the "CoSimulationBaseCouplingSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise NotImplementedError(err_msg)

    ## Check : Called once at the beginning of simulation.
    #          Necessary setup of the solver can be verified here.
    #          IMPORTANT : Check the time step sizes and other setting of
    #                       participating solvers.
    #

    def Check(self):
        for solver_name, solver in self.participating_solvers.items():
           solver.Check()

    ## PrintInfo : Function to display information about the
    #              specifics of the coupled solver.
    #

    def PrintInfo(self):
        cs_tools.PrintInfo("The class", self.__class__.__name__," has the following participants")
        for solver_name, solver in self.participating_solvers.items():
            solver.PrintInfo()

    ## _SynchronizeInputData : Protected Function to obtain new (if any) input data for the
    #                          participating solvers. That will get the data necessary from
    #                          python cosim solver with name solver_name and export (if necessary)
    #                          it to the remote solver of solver_name
    #
    #
    #  So here the _SynchronizeXXXXXData functions exchange data between the python co simulation solvers.
    #  How the solvers get their data from their remote solvers is handled in ImportCouplingInterfaceData and ExportCouplingInterfaceData
    #  functions. Here Import and Export data functions of to_solver and from_solver are called.
    #
    #  @param solver_name     string: name of the solver for which data has to be synchronized
    def _SynchronizeInputData(self, solver_name):
        if self.coupling_started:
            solver = self.participating_solvers[solver_name]
            input_data_list = self.solver_settings[solver_name]["input_data_list"]
            num_input_data = input_data_list.size()

            for i in range(num_input_data):
                input_data = input_data_list[i]
                from_solver = self.participating_solvers[input_data["from_solver"].GetString()]
                data_name = input_data["data_name"].GetString()
                from_solver_data = from_solver.GetInterfaceData(data_name)
                solver_data = solver.GetInterfaceData(input_data["settings"]["destination_data"].GetString())

                if(input_data["settings"].Has("mapper_settings")):
                    solver_data.origin_data = from_solver_data
                    solver_data.mapper_settings = input_data["settings"]["mapper_settings"]

                solver.ImportCouplingInterfaceData(solver_data, from_solver)

    ## _SynchronizeOutputData : Protected Function to synchronize the out put data between the solver
    #                           interface and the remote solver. This assumes that the remote solver
    #                           has output the data in the format specified in the settings (of the out put data def in JSON)
    #
    def _SynchronizeOutputData(self, solver_name):
        if self.coupling_started:
            solver = self.participating_solvers[solver_name]
            output_data_list = self.solver_settings[solver_name]["output_data_list"]
            num_output_data = output_data_list.size()

            for i in range(num_output_data):
                output_data = output_data_list[i]
                to_solver = self.participating_solvers[output_data["to_solver"].GetString()]
                data_name = output_data["data_name"].GetString()
                to_solver_data = to_solver.GetInterfaceData(data_name)
                solver_data = solver.GetInterfaceData(output_data["settings"]["origin_data"].GetString())

                if(output_data["settings"].Has("mapper_settings")):
                    solver_data.destination_data = to_solver_data
                    solver_data.mapper_settings = output_data["settings"]["mapper_settings"]

                solver.ExportCouplingInterfaceData(solver_data, to_solver)

    ## _CreateSolvers : Private Function to make the participating solver objects
    #
    def _CreateSolvers(self, SolversDataMap):
        solvers_map = collections.OrderedDict()
        num_solvers = len(SolversDataMap.keys())
        import KratosMultiphysics.CoSimulationApplication.custom_co_simulation_solver_interfaces.co_simulation_solver_factory as factory

        for solver_name, settings in SolversDataMap.items():
            solver = factory.CreateSolverInterface(solver_name,cs_data_structure.Parameters(settings))
            solvers_map[solver_name] = solver

        return solvers_map

    ## _GetSolverCoSimulationDetails : Private Function to obtain a dict of setting with solver
    #                                  name as key
    #
    def _GetSolverCoSimulationDetails(self,co_simulation_solver_settings):
        num_solvers = co_simulation_solver_settings.size()
        solver_cosim_details = {}
        for i_solver in range(num_solvers):
            solver_name = co_simulation_solver_settings[i_solver]["name"].GetString()
            solver_cosim_details[solver_name] = co_simulation_solver_settings[i_solver]
        # TODO check if the data is consistently defined! => maybe do at another place though...
        # - input in one is output in another
        # - one IO is defined for each data_name
        # - if the same data is defined multiple times
        # - check if data format has been specified
        return solver_cosim_details

    ## _CreateFilters : Protected Function to make filter objects list and store in the datafield
    #
    #  @param conv_acc_settings dict: setting of the convergence accelerator to be make
    def _CreateFilters(self, co_simulation_solver_settings): # probably better in some utils file
        import KratosMultiphysics.CoSimulationApplication.custom_convergence_accelerators.co_simulation_convergence_accelerator_factory as factory
        num_solvers = co_simulation_solver_settings.size()
        solver_cosim_details = {}
        for i_solver in range(num_solvers):
            solver_name = co_simulation_solver_settings[i_solver]["name"].GetString()
            solver = self.participating_solvers[solver_name]

            ## for all the input data
            input_data_list = self.solver_settings[solver_name]["input_data_list"]
            num_input_data = input_data_list.size()
            for i in range(num_input_data):
                input_data = input_data_list[i]
                filters_list = input_data["filters"]
                destination_data = solver.GetInterfaceData(input_data["settings"]["destination_data"].GetString())
                for filter in filters_list:
                    accelerator = factory.CreateConvergenceAccelerator(accelerator_settings, solver) ## TODO: should change to filter
                    destination_data.filters.append(accelerator)


            ## for all the output data
            output_data_list = self.solver_settings[solver_name]["output_data_list"]
            num_output_data = output_data_list.size()
            for i in range(num_output_data):
                output_data = output_data_list[i]
                filters_list = input_data["filters"]
                origin_data = solver.GetInterfaceData(input_data["settings"]["origin_data"].GetString())
                for filter in filters_list:
                    accelerator = factory.CreateConvergenceAccelerator(accelerator_settings, solver) ## TODO: should change to filter
                    origin_data.filters.append(accelerator)

    ## _CreateConvergenceCriteria : Private Function to make convergence criteria objects list #Comment protected
    #
    #  @param conv_acc_settings dict: setting of the convergence criteria to be make
    def _CreateConvergenceCriteria(self, conv_criteria_settings): # probably better in some utils file
        conv_criteria = []
        import KratosMultiphysics.CoSimulationApplication.base_co_simulation_classes.co_simulation_base_convergence_criteria as criteria
        num_criteria = conv_criteria_settings.size()
        for i in range(num_criteria):
            criteria_setting = conv_criteria_settings[i]
            solver_name = criteria_setting["solver"].GetString()
            solver = self.participating_solvers[solver_name]
            criteria = criteria.Create(criteria_setting, solver) # Change to use interface data
            conv_criteria.append(criteria)

        return conv_criteria