# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.factories.io_factory as io_factory
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities

def Create(settings, name):
    raise Exception('"CoSimulationSolverWrapper" is a baseclass and cannot be used directly!')

class CoSimulationSolverWrapper:
    """Baseclass for the solver wrappers used for CoSimulation
    It wraps solvers used in the CoSimulation
    """
    def __init__(self, settings, model, solver_name):
        """Constructor of the Base Solver Wrapper

        The derived classes should do the following things in their constructors:
        1. call the base-class constructor (i.e. the constructor of this class => CoSimulationSolverWrapper)
        2. create the ModelParts required for the CoSimulation
        3. Optional: call "_AllocateHistoricalVariablesFromCouplingData" to allocate the nodal historical variables on the previously created ModelParts (this should not be necessary for Kratos)
           => this has to be done before the meshes/coupling-interfaces are read/received/imported (due to how the memory allocation of Kratos works for historical nodal values)
        """

        # Every SolverWrapper has its own model, because:
        # - the names can be easily overlapping (e.g. "Structure.Interface")
        # - Solvers should not be able to access the data of other solvers directly!
        self.model = model
        if self.model == None:
            self.model = KM.Model()
        elif not isinstance(self.model, KM.Model):
            err_msg  = 'A solver wrapper can either be passed a Model\n'
            err_msg += 'or None, got object of type "{}"'.format(type(self.model))
            raise Exception(err_msg)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultParameters())

        self.name = solver_name
        if "." in self.name:
            raise Exception("cannot contain dot!")

        self.echo_level = self.settings["echo_level"].GetInt()

        self.data_communicator = self._GetDataCommunicator()

        # The IO is only used if the corresponding solver is used in coupling and it initialized from the "higher instance, i.e. the coupling-solver
        self.__io = None

    def _GetSolver(self, solver_name):
        raise Exception('Trying to get SolverWrapper "{}" of "{}" which is not a coupled solver!'.format(solver_name, self.name))

    def Initialize(self):
        self.data_dict = {data_name : CouplingInterfaceData(data_config, self.model, data_name, self.name) for (data_name, data_config) in self.settings["data"].items()}

        if self.__HasIO():
            self.__GetIO().Initialize()

    def Finalize(self):
        if self.__HasIO():
            self.__GetIO().Finalize()

    def AdvanceInTime(self, current_time):
        # in case a solver does not provide time information (e.g. external or steady solvers),
        # then this solver should return "0.0" here
        raise Exception('"AdvanceInTime" must be implemented in the derived class!')

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


    def CreateIO(self, io_echo_level):
        if self.__HasIO():
            raise Exception('IO for solver "{}" is already created!'.format(self.name))

        io_settings = self.settings["io_settings"]

        if not io_settings.Has("echo_level"):
            io_settings.AddEmptyValue("echo_level").SetInt(io_echo_level)

        self.__io = io_factory.CreateIO(self.settings["io_settings"], self.model, self.name, self.data_communicator, self._GetIOType())

    def ImportCouplingInterface(self, interface_config):
        if self.echo_level > 2:
            cs_tools.cs_print_info("CoSimulationSolverWrapper", 'Importing coupling interface "{}" of solver: "{}"'.format(colors.magenta(interface_config["model_part_name"]), colors.blue(self.name)))
        self.__GetIO().ImportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        if self.echo_level > 2:
            cs_tools.cs_print_info("CoSimulationSolverWrapper", 'Exporting coupling interface "{}" of solver: "{}"'.format(colors.magenta(interface_config["model_part_name"]), colors.blue(self.name)))
        self.__GetIO().ExportCouplingInterface(interface_config)

    def ImportData(self, data_config):
        if self.echo_level > 2:
            cs_tools.cs_print_info("CoSimulationSolverWrapper", 'Importing data of solver: "{}" with type: "{}"'.format(colors.blue(self.name), data_config["type"]))
        self.__GetIO().ImportData(data_config)

    def ExportData(self, data_config):
        if self.echo_level > 2:
            cs_tools.cs_print_info("CoSimulationSolverWrapper", 'Exporting data of solver: "{}" with type: "{}"'.format(colors.blue(self.name), data_config["type"]))
        self.__GetIO().ExportData(data_config)

    def IsDefinedOnThisRank(self):
        return self.data_communicator.IsDefinedOnThisRank()

    def GetInterfaceData(self, data_name):
        try:
            return self.data_dict[data_name]
        except KeyError:
            raise Exception('Requested data field "{}" does not exist for solver "{}"'.format(data_name, self.name))

    def PrintInfo(self):
        if self.echo_level > 0 and self.data_communicator.IsDefinedOnThisRank():
            cs_tools.cs_print_info("CoSimulationSolverWrapper", self._ClassName(), self.name, "Using {} mpi-processes".format(self.data_communicator.Size()))

    def Check(self):
        cs_tools.cs_print_warning("CoSimulationSolverWrapper", "your solver does not implement Check!!!")

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    def _GetIOType(self):
        # only external solvers have to specify sth here / override this
        return "dummy_io"

    def __GetIO(self):
        if not self.__HasIO():
            raise Exception('IO for solver "{}" is not created!'.format(self.name))
        return self.__io

    def __HasIO(self):
        return self.__io is not None

    def _GetDataCommunicator(self):
        if len(self.settings["mpi_settings"].keys()) > 0:
            if self.settings["mpi_settings"]["num_processes"].GetInt() == 1:
                return data_communicator_utilities.GetRankZeroDataCommunicator()
            return data_communicator_utilities.CreateDataCommunicatorWithNProcesses(self.settings["mpi_settings"])
        else:
            # if no special input is specified use the default implementation from the baseclass
            return KM.ParallelEnvironment.GetDefaultDataCommunicator()

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "type"                    : "",
            "solver_wrapper_settings" : {},
            "io_settings"             : {},
            "data"                    : {},
            "mpi_settings"            : {},
            "echo_level"              : 0
        }""")
