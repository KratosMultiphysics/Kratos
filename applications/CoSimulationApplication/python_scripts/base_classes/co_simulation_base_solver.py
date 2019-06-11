from __future__ import print_function, absolute_import, division

# Other imports
import KratosMultiphysics.CoSimulationApplication.solver_interfaces.co_simulation_io_factory as io_factory

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
# Other imports
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

def CreateSolver(model, settings):
    return CoSimulationBaseSolver(model, settings)

class CoSimulationBaseSolver(object):
    """The base class for the CoSimulation Solvers
    """
    def __init__(self, model, settings, solver_name):
        """Constructor of the Base-Solver
        Deriving classes should call it in their constructors
        """

        self.model = model

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())

        self.name = solver_name
        self.echo_level = self.settings["echo_level"].GetInt()
        self.io_is_initialized = False
        self.data_dict = self.__CreateInterfaceDataDict()

    def Initialize(self):
        pass

    def InitializeIO(self, solvers, io_echo_level):
        if self.io_is_initialized:
            raise Exception('IO for "' + self.name + '" is already initialized!')

        self.io = io_factory.CreateIO(self._GetIOName(),
                                      self.model,
                                      self.settings["io_settings"])
        self.io.SetEchoLevel(io_echo_level)
        self.io_is_initialized = True

    def Finalize(self):
        pass

    def AdvanceInTime(self, current_time):
        return current_time + self.settings["time_step"] # needed if this solver is used as dummy

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

    def GetInterfaceData(self, data_name):
        try:
            return self.data_dict[data_name]
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

    ## __CreateInterfaceDataDict : Private Function to obtain the map of data objects
    #
    #  @param self            The object pointer.
    def __CreateInterfaceDataDict(self):
        data_dict = dict()
        for data_name, data_config in self.settings["data"].items():
            data_dict[data_name] = CouplingInterfaceData(data_config, self.model)

        return data_dict

    @classmethod
    def _GetDefaultSettings(cls):
        return cs_tools.cs_data_structure.Parameters("""{
            "type"        : "",
            "io_settings" : {},
            "settings"    : {},
            "data"        : {},
            "echo_level"  : 0
        }""")
