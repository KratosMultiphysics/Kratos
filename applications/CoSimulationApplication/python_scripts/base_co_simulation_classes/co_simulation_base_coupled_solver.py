# co simulation imports
import co_simulation_tools as tools
# Importing the CoSimulation application
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
from co_simulation_base_solver import CoSimulationBaseSolver

# Other imports


##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
#  This class is intended to server as the base class for all the coupled solvers.
class CoSimulationBaseCoupledSolver(CoSimulationBaseSolver):
    ## The constructor
    #
    #  @param self            The object pointer.
    #  @param custom_settings     parameters for configuring the CoSimulationBaseCoupledSolver
    def __init__(self, custom_settings):

        ##settings string in json format
        # default values for all available settings
        # for mandatory settings, the type is defined
        self.full_settings = custom_settings
        #super(CoSimulationBaseCoupledSolver,self).__init__(custom_settings)
        defaultSettings = {}
        defaultSettings["echo_level"] = 1
        defaultSettings["max_coupling_iterations"] = 10
        defaultSettings["participants"] = list
        defaultSettings["start_coupling_time"] = float #MANDATORY
        self.settings = tools.ValidateAndAssignInputParameters(defaultSettings, custom_settings["coupled_solver_settings"], False)
        self.number_of_participants = len( self.settings['participants'] )
        self.echo_level = self.settings["echo_level"]

        # Get the participating solvers a map with their names and objects
        self.participating_solvers = self._GetSolvers(self.full_settings['solvers'])

        self.solver_settings = self._GetSolverCoSimulationDetails(self.full_settings['coupled_solver_settings']["participants"])

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if "start_coupling_time" in self.settings:
            self.start_coupling_time = self.settings["start_coupling_time"]
        if self.start_coupling_time > 0.0:
            self.coupling_started = False
        else:
            self.coupling_started = True

    ## Initialize : Initialize function. Called only once
    #               all member variables are initialized here.
    #  @param self            The object pointer.
    def Initialize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Initialize()

    ## Finalize : Finalize function. Called only once
    #               all member variables are finalized here.
    #  @param self            The object pointer.
    def Finalize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Finalize()

    ## Predict : Predict the solution of the next solution step.
    #
    #  @param self            The object pointer.
    def Predict(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Predict()

    ## InitializeSolutionStep : Called once in the beginning of the solution step
    #
    #  @param self            The object pointer.
    def InitializeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.InitializeSolutionStep()

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #
    #  @param self            The object pointer.
    def FinalizeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.FinalizeSolutionStep()

    ## OutputSolutionStep : Called once at the end of the solution step.
    #                       The output of the solvers and / or cosimulation output
    #                       can be performed in this function.
    #
    #  @param self            The object pointer.
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
    #  @param self            The object pointer.
    def SolveSolutionStep(self):
        err_msg  = 'Calling "SolveSolutionStep" of the "CoSimulationBaseCouplingSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise Exception(err_msg)

    ## Check : Called once at the beginning of simulation.
    #          Necessary setup of the solver can be verified here.
    #          IMPORTANT : Check the time step sizes and other setting of
    #                       participating solvers.
    #
    #  @param self            The object pointer.
    def Check(self):
        for solver_name, solver in self.participating_solvers.items():
           solver.Check()

    ## PrintInfo : Function to display information about the
    #              specifics of the coupled solver.
    #
    #  @param self            The object pointer.
    def PrintInfo(self):
        print("The class", self.__class__.__name__," has the following participants:")
        for solver_name, solver in self.participating_solvers.items():
            solver.PrintInfo()

    ## _SynchronizeInputData : Private Function to obtain new (if any) input data for the
    #                          participating solvers. That will get the data necessary from
    #                          python cosimm solver with name solver_name and export (if necessary)
    #                          it to the remote solver of solver_name
    #
    #
    #  So here the _SynchronizeXXXXXData functions exchange data between the python co simulation solvers.
    #  How the solvers get their data from their remote solvers is handled in ImportData and ExportData
    #  functions. Here Import and Export data functions of to_solver and from_solver are called.
    #
    #  @param self            The object pointer.
    #  @param solver_name     string: name of the solver for which data has to be synchronized
    def _SynchronizeInputData(self, solver_name):
        if self.coupling_started:
            solver = self.participating_solvers[solver_name]
            input_data_list = self.solver_settings[solver_name]["input_data_list"]
            for input_data in input_data_list:
                from_solver = self.participating_solvers[input_data["from_solver"]]
                data_name = input_data["data_name"]
                from_solver_data_conf = from_solver.GetDataConfig(data_name)
                current_sovler_data_conf = solver.GetDataConfig(data_name)
                if( from_solver_data_conf["dimension"] == current_sovler_data_conf["dimension"] ):
                    solver.ImportData(data_name, from_solver)

    ## _SynchronizeOutputData : Private Function to synchronize the out put data between the solver
    #                           interface and the remote solver. This assumes that the remote solver
    #                           has output the data inthe format specified in the settings (of the out put data def in JSON)
    #
    #  @param self            The object pointer.
    def _SynchronizeOutputData(self, solver_name):
        if self.coupling_started:
            solver = self.participating_solvers[solver_name]
            output_data_list = self.solver_settings[solver_name]["output_data_list"]
            for output_data in output_data_list:
                data_name = output_data["data_name"]
                data_conf = solver.GetDataConfig(data_name)
                solver.ExportData(data_conf)

    ## _GetSolvers : Private Function to make the paticipating solver objects
    #
    #  @param self            The object pointer.
    def _GetSolvers(self, SolversDataList):
        solvers_map = {}
        num_solvers = len(SolversDataList)

        for i in range(0,num_solvers):
            import co_simulation_solver_factory as factory
            solver = factory.CreateSolverInterface(SolversDataList[i])
            solver_name = SolversDataList[i]["name"]
            solvers_map[solver_name] = solver

        return solvers_map

    ## _GetSolverCoSimulationDetails : Private Function to obtain a dict of setting with solver
    #                                  name as key
    #
    #  @param self            The object pointer.
    def _GetSolverCoSimulationDetails(self,co_simulation_solver_settings):
        num_solvers = len(co_simulation_solver_settings)
        solver_cosim_details = {}
        for solver_settings in co_simulation_solver_settings:
            solver_name = solver_settings["name"]
            solver_cosim_details[solver_name] = solver_settings
        # TODO check if the data is consitently defined! => maybe do at another place though...
        # - input in one is output in another
        # - one IO is defined for each data_name
        # - if the same data is defined multiple times
        # - check if data format has been specified
        return solver_cosim_details