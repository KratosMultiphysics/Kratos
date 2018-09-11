from __future__ import print_function, absolute_import, division

# Other imports
from . import io_factory

def CreateSolver(cosim_solver_settings, level):
    return CoSimulationBaseSolver(cosim_solver_settings, level)

class CoSimulationBaseSolver(object):
    ## Constructor :  The base class for the CoSimulation Solver interfaces
    #                  Constructor of the Base-Solver interface
    #                  Deriving classes should call it in their constructors
    #
    #  @param self                      The object pointer.
    #  @param cosim_solver_settings     python dictionary : with the solver settings.
    def __init__(self, cosim_solver_settings):
        self.cosim_solver_settings = cosim_solver_settings
        self.echo_level = 0
        if "echo_level" in self.cosim_solver_settings:
            self.echo_level = self.cosim_solver_settings["echo_level"]
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
    def InitializeIO(self, io_echo_level):
        solver_name = self.cosim_solver_settings["name"]
        if self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is already initialized!')

        self.io = io_factory.CreateIO(self._GetIOName(),
                                      solvers,
                                      solver_name)
        self.io.SetEchoLevel(io_echo_level)
        self.io_is_initialized = True

    ## Finalize :  Initialize function for the solver class.
    #               finalization of the variables and objects to be done here.
    #  @param self                      The object pointer.
    def Finalize(self):
        pass

    def AdvanceInTime(self, current_time):
        pass

    def Predict(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def OutputSolutionStep(self):
        pass

    def SolveSolutionStep(self):
        pass

    def ImportData(self, data_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportData(data_name, from_client)
    def ImportMesh(self, mesh_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ImportMesh(mesh_name, from_client)

    def ExportData(self, data_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportData(data_name, to_client)
    def ExportMesh(self, mesh_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + solver_name + '" is not initialized!')
        self.io.ExportMesh(mesh_name, to_client)

    def GetDataDefinition(self, data_name):
        return self.cosim_solver_settings["data"][data_name]

    def GetDeltaTime(self):
        raise Exception('!!!WARNING!!!  Calling "GetDeltaTime" function from base Co-Simulation Solver class!')

    def PrintInfo(self):
        '''
        This function can be filled if desired, e.g. to print settings at higher echo-levels
        '''
        pass

    def SetEchoLevel(self, level):
        self.echo_level = level

    def Check(self):
        print("!!!WARNING!!! Calling Check from base Co-Simulation Solver class !!!")

    def _GetIOType(self):
        raise Exception('"_GetIOName" function must be implemented in derived class!')

