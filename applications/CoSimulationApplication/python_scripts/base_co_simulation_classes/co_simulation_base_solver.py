from __future__ import print_function, absolute_import, division

# Other imports
import co_simulation_io_factory as io_factory
import co_simulation_tools as tools

##
#  IMPORTANT : This is a BASE CLASS
#               Please do not change any thing in this class.
#
# This Class servers as a base class for all the Solver interfaces which will be implemented
class CoSimulationBaseSolver(object):
    ## Constructor :  The base class for the CoSimulation Solver interfaces
    #                  Constructor of the Base-Solver interface
    #                  Deriving classes should call it in their constructors
    #
    #  @param self                      The object pointer.
    #  @param cosim_solver_settings     python dictionary : with the solver settings.
    def __init__(self, cosim_solver_settings):
        defaultSettings = {}
        defaultSettings["name"] = str # MANDATORY
        defaultSettings["settings"] = dict # MANDATORY
        defaultSettings["data"] = list # MANDATORY
        defaultSettings["echo_level"] = 0
        self.cosim_solver_settings = tools.ValidateAndAssignInputParameters(defaultSettings, cosim_solver_settings, False)
        self.SetEchoLevel( self.cosim_solver_settings["echo_level"] )

        self.data_list = self._GetDataList()

        self.io_is_initialized = False

    ## Initialize : Initialize function for the solver class. Necessary
    #               initialization of the variables and objects to be done here.
    #  @param self                      The object pointer.
    def Initialize(self):
        pass

    ## InitializeIO : Initialize the IO class for this solver.
    #                   usually a particular type of solver has a particular default IO type
    #  @param self                      The object pointer.
    #  @param io_echo_level             int : echo level for the io to be initialized.
    def InitializeIO(self, io_echo_level=0):
        solver_name = self.cosim_solver_settings["name"]
        if self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is already initialized!')

        self.io = io_factory.CreateIO(self._GetIOName(), solver_name)
        self.io.SetEchoLevel(io_echo_level)
        self.io_is_initialized = True

    ## Finalize :  Initialize function for the solver class.
    #               finalization of the variables and objects to be done here.
    #  @param self                      The object pointer.
    def Finalize(self):
        pass

    ## AdvanceInTime :  This function will advance the current solver to the given current_time
    #
    #  @param self                      The object pointer.
    #  @current_time                    The time to which the current solver should advance
    def AdvanceInTime(self, current_time):
        pass

    ## Predict : Predict the solution of the next solution step.
    #
    #  @param self            The object pointer.
    def Predict(self):
        pass

    ## InitializeSolutionStep : Called once in the beginning of the solution step
    #
    #  @param self            The object pointer.
    def InitializeSolutionStep(self):
        pass

    ## FinalizeSolutionStep : Called once at the end of the solution step
    #
    #  @param self            The object pointer.
    def FinalizeSolutionStep(self):
        pass

    ## OutputSolutionStep : Called once at the end of the solution step.
    #                       The output of the solvers and / or cosimulation output
    #                       can be performed in this function.
    #
    #  @param self            The object pointer.
    def OutputSolutionStep(self):
        pass

    ## SolveSolutionStep : This function implements the coupling workflow
    #
    #  @param self            The object pointer.
    def SolveSolutionStep(self):
        pass

    ## ImportData : This function imports the requesed data from
    #               from_client
    #
    #  @param self            The object pointer.
    #  @param data_name       string : Name of the data to be imported from from_client
    #  @param from_client     python obj : The client from which data_name has to be imported
    def ImportData(self, data_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        data_conf = self.GetDataConfig(data_name)
        self.io.ImportData(data_conf, from_client)

    ## ImportMesh : This function imports the requesed surface/volume
    #               mesh from from_client
    #
    #  @param self            The object pointer.
    #  @param mesh_name       string : Name of the mesh to be imported from from_client
    #  @param from_client     python obj : The client from which data_name has to be imported
    def ImportMesh(self, mesh_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportMesh(mesh_name, from_client)

    ## ExportData : This function exports the requesed data to
    #               to_client
    #
    #  @param self            The object pointer.
    def ExportData(self, data_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportData(data_name, to_client)

    ## ExportMesh : This function exports the requesed surface/volume
    #               to to_client
    #
    #  @param self            The object pointer.
    def ExportMesh(self, mesh_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportMesh(mesh_name, to_client)

    ## GetDataDefinition : Function to get the data configuration belonging
    #                      requested data field
    #
    #  @param self            The object pointer.
    def GetDataDefinition(self, data_name):
        return self.cosim_solver_settings["data"][data_name]

    ## GetDeltaTime : Function to obtain the time step of this solver
    #
    #  @param self            The object pointer.
    def GetDeltaTime(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling GetDeltaTime function from base Co-Simulation Solver class!" + tools.bcolors.ENDC)

    ## PrintInfo : Function to display information about the
    #              specifics of the coupled solver.
    #
    #  @param self            The object pointer.
    def PrintInfo(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling PrintInfo function from base Co-Simulation Solver class!" + tools.bcolors.ENDC)

    ## SetEchoLevel : Function to set the Echo level for this solver
    #
    #  @param self            The object pointer.
    def SetEchoLevel(self, level):
        self.echo_level = level

    ## Check : Called once at the beginning of simulation.
    #          Necessary setup of the solver can be verified here.
    #
    #  @param self            The object pointer.
    def Check(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling Check function from base Co-Simulation Solver class!" + tools.bcolors.ENDC)

    ## _GetIOType : Private Function to obtain the type of IO this solver has as default
    #
    #  @param self            The object pointer.
    def _GetIOType(self):
        raise NotImplementedError(tools.bcolors.FAIL + "CoSimulationBaseSolver : Calling _GetIOName function from base Co-Simulation Solver class!" + tools.bcolors.ENDC)

    ## _GetIOType : Private Function to obtain the list of data objects
    #
    #  @param self            The object pointer.
    def _GetDataList(self):
        data_list = {}
        for data_def in self.cosim_solver_settings["data"]:
            data_conf = tools.GetDataConfig(data_def)
            data_list[data_conf["name"]] = data_conf

        return data_list
