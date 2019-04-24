from __future__ import print_function, absolute_import, division

# Other imports
import KratosMultiphysics.CoSimulationApplication.solver_interfaces.co_simulation_io_factory as io_factory

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as tools
# Other imports
cs_data_structure = tools.cs_data_structure
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

def CreateSolver(model, cosim_solver_settings):
    return CoSimulationBaseSolver(model, cosim_solver_settings)

class CoSimulationBaseSolver(object):
    """The base class for the CoSimulation Solvers
    """
    def __init__(self, model, cosim_solver_settings, solver_name):
        """Constructor of the Base-Solver
        Deriving classes should call it in their constructors
        """
        default_settings = cs_data_structure.Parameters("""{
            "solver_type" : "",
            "io_settings" : {},
            "settings"    : {},
            "data"        : [],
            "echo_level"  : 0
        }""")

        self.model = model
        self.cosim_solver_settings = cosim_solver_settings
        self.cosim_solver_settings.AddMissingParameters(default_settings)
        self.name = solver_name
        self.echo_level = self.cosim_solver_settings["echo_level"].GetInt()
        self.io_is_initialized = False
        self.data_map = self.__CreateInterfaceDataMap()

    def Initialize(self):
        pass

    def InitializeIO(self, solvers, io_echo_level):
        if self.io_is_initialized:
            raise Exception('IO for "' + self.name + '" is already initialized!')

        self.io = io_factory.CreateIO(self._GetIOName(),
                                      self.model,
                                      self.cosim_solver_settings["io_settings"])
        self.io.SetEchoLevel(io_echo_level)
        self.io_is_initialized = True

    def Finalize(self):
        pass

    def AdvanceInTime(self, current_time):
        return current_time + self.cosim_solver_settings["time_step"] # needed if this solver is used as dummy

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

    def ImportCouplingInterfaceData(self, data_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + self.name + '" is not initialized!')
        self.io.ImportCouplingInterfaceData(data_name, from_client)
    def ImportCouplingInterface(self, geometry_name, from_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + self.name + '" is not initialized!')
        self.io.ImportCouplingInterface(geometry_name, from_client)

    def ExportCouplingInterfaceData(self, data_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + self.name + '" is not initialized!')
        self.io.ExportCouplingInterfaceData(data_name, to_client)
    def ExportCouplingInterface(self, geometry_name, to_client):
        if not self.io_is_initialized:
            raise Exception('IO for "' + self.name + '" is not initialized!')
        self.io.ExportCouplingInterface(geometry_name, to_client)

    def GetInterfaceDataOld(self, data_name):
        return self.cosim_solver_settings["data"][data_name]


    def GetInterfaceData(self, data_name):
        try:
            return self.data_map[data_name]
        except KeyError:
            raise Exception("Requested data field " + data_name + " does not exist in the solver ")

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

    def IsDistributed(self):
        '''Returns whether this solver is executed distributed Aka MPI-parallel
        '''
        return False

    def _GetIOName(self):
        raise Exception('"_GetIOName" function must be implemented in derived class!')

    ## __CreateInterfaceDataMap : Private Function to obtain the map of data objects
    #
    #  @param self            The object pointer.
    def __CreateInterfaceDataMap(self):
        data_map = dict()
        # # num_data = self.cosim_solver_settings["data"].size()
        # print(self.cosim_solver_settings["data"].keys())
        # errrr

        print("\n\nsolver_name", self.name)
        # print(self.cosim_solver_settings.PrettyPrintJsonString())
        # for i in range(self.cosim_solver_settings["data"].size()):
        for i in range(self.cosim_solver_settings["data"].size()):
            data_conf = self.cosim_solver_settings["data"][i]
            data_name = data_conf["name"].GetString()
            print(data_conf.PrettyPrintJsonString())
            print(data_name)
            data_obj = CouplingInterfaceData(data_conf, self)
            data_map[data_name] = data_obj

        return data_map
