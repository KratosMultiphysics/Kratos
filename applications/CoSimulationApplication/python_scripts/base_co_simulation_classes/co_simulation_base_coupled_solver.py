# co simulation imports
import co_simulation_tools as tools
from base_co_simulation_classes.co_simulation_base_solver import CoSimulationBaseSolver
# Other imports
import co_simulation_data_structure
cs_data_structure = co_simulation_data_structure.__DATA_STRUCTURE__
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
        self.settings.ValidateAndAssignDefaults(default_setting)
        self.number_of_participants = self.settings['participants'].size()
        self.echo_level = self.settings["echo_level"].GetInt()
        self.num_coupling_iterations = self.settings["num_coupling_iterations"].GetInt()

        # Get the participating solvers a map with their names and objects
        self.participating_solvers = self._GetSolvers(self.full_settings['solvers'])

        self.solver_settings = self._GetSolverCoSimulationDetails(self.full_settings['coupled_solver_settings']["participants"])

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
        raise Exception(err_msg)

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
        print("The class", self.__class__.__name__," has the following participants:")
        for solver_name, solver in self.participating_solvers.items():
            solver.PrintInfo()

    ## _SynchronizeInputData : Private Function to obtain new (if any) input data for the
    #                          participating solvers. That will get the data necessary from
    #                          python cosim solver with name solver_name and export (if necessary)
    #                          it to the remote solver of solver_name
    #
    #
    #  So here the _SynchronizeXXXXXData functions exchange data between the python co simulation solvers.
    #  How the solvers get their data from their remote solvers is handled in ImportData and ExportData
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
                from_solver_data_conf = from_solver.GetDataConfig(data_name)
                solver_data_conf = solver.GetDataConfig(input_data["settings"]["destination_data"].GetString())
                if(input_data["settings"].Has("mapper_settings")):
                    solver_data_conf.AddValue("origin_data_config",from_solver_data_conf)
                    solver_data_conf.AddValue("mapper_settings", input_data["settings"]["mapper_settings"])

                solver.ImportData(solver_data_conf, from_solver)

    ## _SynchronizeOutputData : Private Function to synchronize the out put data between the solver
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
                data_name = output_data["data_name"].GetString()
                data_conf = solver.GetDataConfig(data_name)
                solver.ExportData(data_conf)

    ## _GetSolvers : Private Function to make the participating solver objects
    #
    def _GetSolvers(self, SolversDataMap):
        solvers_map = collections.OrderedDict()
        num_solvers = len(SolversDataMap.keys())
        import custom_co_simulation_solver_interfaces.co_simulation_solver_factory as factory

        for solver_name, settings in SolversDataMap.items():
            solver = factory.CreateSolverInterface(solver_name,settings)
            solvers_map[solver_name] = solver

        return solvers_map

    ## _GetSolverCoSimulationDetails : Private Function to obtain a dict of setting with solver
    #                                  name as key
    #
    def _GetSolverCoSimulationDetails(self,co_simulation_solver_settings):
        num_solvers = co_simulation_solver_settings.size()
        solver_cosim_details = {}

        for i in range(num_solvers):
            solver_settings = co_simulation_solver_settings[i]
            solver_name = solver_settings["name"].GetString()
            solver_cosim_details[solver_name] = solver_settings
        # TODO check if the data is consistently defined! => maybe do at another place though...
        # - input in one is output in another
        # - one IO is defined for each data_name
        # - if the same data is defined multiple times
        # - check if data format has been specified
        return solver_cosim_details

    ## _GetConvergenceAccelerators : Private Function to make convergence accelerator objects list
    #
    #  @param conv_acc_settings dict: setting of the convergence accelerator to be make
    def _GetConvergenceAccelerators(self, conv_acc_settings):
        conv_accelerators = []
        num_acceleratos = conv_acc_settings.size()
        import custom_convergence_accelerators.co_simulation_convergence_accelerator_factory as factory
        for i in range(num_acceleratos):
            accelerator_settings = conv_acc_settings[i]
            solver_name = accelerator_settings["data"]["solver"].GetString()
            solver = self.participating_solvers[solver_name]
            accelerator = factory.CreateConvergenceAccelerator(accelerator_settings, solver)
            conv_accelerators.append(accelerator)

        return conv_accelerators

    ## _GetConvergenceCriteria : Private Function to make convergence criteria objects list
    #
    #  @param conv_acc_settings dict: setting of the convergence criteria to be make
    def _GetConvergenceCriteria(self, conv_criteria_settings):
        conv_criteria = []
        import base_co_simulation_classes.co_simulation_base_convergence_criteria as criteria
        num_criteria = conv_criteria_settings.size()
        for i in range(num_criteria):
            criteria_setting = conv_criteria_settings[i]
            solver_name = criteria_setting["solver"].GetString()
            solver = self.participating_solvers[solver_name]
            criteria = criteria.Create(criteria_setting, solver)
            conv_criteria.append(criteria)

        return conv_criteria