from __future__ import print_function, absolute_import, division

# Other imports
import co_simulation_ios.co_simulation_io_factory as io_factory

def CreateSolver(cosim_solver_settings, level):
    return CoSimulationBaseSolver(cosim_solver_settings, level)

class CoSimulationBasePhysicsSolver(object):
    """The base class for the CoSimulation Solvers
    The intention is that every solver that derives from this class
    can be used standalone.
    """
    def __init__(self, cosim_solver_settings, level):
        """Constructor of the Base-Solver
        Deriving classes should call it in their constructors
        """
        self.cosim_solver_settings = cosim_solver_settings
        self.lvl = level
        self.echo_level = 0
        if "echo_level" in self.cosim_solver_settings:
            self.echo_level = self.cosim_solver_settings["echo_level"]
        self.io_is_initialized = False

    def Initialize(self):
        if self.__IsSolvingRank():
            self._Initialize()

    def InitializeIO(self, solvers, io_echo_level):
        if self.__IsSolvingRank():
            solver_name = self.cosim_solver_settings["name"]
            if self.io_is_initialized:
                raise Exception('IO for "' + solver_name + '" is already initialized!')

            self.io = io_factory.CreateIO(self._GetIOName(),
                                        solvers,
                                        solver_name,
                                        self.lvl)
            self.io.SetEchoLevel(io_echo_level)
            self.io_is_initialized = True

    def Finalize(self):
        if self.__IsSolvingRank():
            self._Finalize()

    def AdvanceInTime(self, current_time):
        return current_time + self.cosim_solver_settings["time_step"] # needed if this solver is used as dummy

    def Predict(self):
        if self.__IsSolvingRank():
            self._Predict()

    def InitializeSolutionStep(self):
        if self.__IsSolvingRank():
            self._InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        if self.__IsSolvingRank():
            self._FinalizeSolutionStep()

    def OutputSolutionStep(self):
        if self.__IsSolvingRank():
            self._OutputSolutionStep()

    def SolveSolutionStep(self):
        if self.__IsSolvingRank():
            self._SolveSolutionStep()

    def ImportData(self, data_name, from_client):
        if self.__IsSolvingRank():
            if not self.io_is_initialized:
                raise Exception('IO for "' + solver_name + '" is not initialized!')
            self.io.ImportData(data_name, from_client)
    def ImportMesh(self, mesh_name, from_client):
        if self.__IsSolvingRank():
            if not self.io_is_initialized:
                raise Exception('IO for "' + solver_name + '" is not initialized!')
            self.io.ImportMesh(mesh_name, from_client)

    def ExportData(self, data_name, to_client):
        if self.__IsSolvingRank():
            if not self.io_is_initialized:
                raise Exception('IO for "' + solver_name + '" is not initialized!')
            self.io.ExportData(data_name, to_client)
    def ExportMesh(self, mesh_name, to_client):
        if self.__IsSolvingRank():
            if not self.io_is_initialized:
                raise Exception('IO for "' + solver_name + '" is not initialized!')
            self.io.ExportMesh(mesh_name, to_client)

    def GetDataDefinition(self, data_name):
        return self.cosim_solver_settings["data"][data_name]

    def GetBufferSize(self):
        raise Exception('"GetBufferSize" function must be implemented in derived class!')

    def GetDeltaTime(self):
        raise Exception('"GetDeltaTime" function must be implemented in derived class!')

    def PrintInfo(self):
        '''This function can be filled if desired, e.g. to print settings at higher echo-levels
        '''
        pass

    def SetEchoLevel(self, level):
        self.echo_level = level

    def Check(self):
        print("!!!WARNING!!! your solver does not implement Check!!!")

    def _GetIOName(self):
        raise Exception('"_GetIOName" function must be implemented in derived class!')

    def __IsSolvingRank(self):
        return self.IsDistributed or (co_simulation_tools.COSIM_SPACE.Rank() == 0)