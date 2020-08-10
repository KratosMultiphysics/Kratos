# Importing the Kratos Library
import KratosMultiphysics as KM

def Create(settings, model, solver_name):
    raise Exception('"CoSimulationIO" is a baseclass and cannot be used directly!')

class CoSimulationIO(object):
    """Baseclass defining the interface for the input and output methods
    for the communication with external solvers
    """
    def __init__(self, settings, model, solver_name):
        self.model = model
        self.solver_name = solver_name # name of the owning solver

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def ImportCouplingInterface(self, interface_config):
        """Imports coupling interface from an external solver
        External solver sends, CoSimulation receives

        @param interface_config <python dictionary> : configuration of the interface to be imported
        """
        raise NotImplementedError("This function has to be implemented in the derived class!")

    def ExportCouplingInterface(self, interface_config):
        """Exports coupling interface to an external solver
        CoSimulation sends, external solver receives

        @param interface_config <python dictionary> : configuration of the interface to be exported
        """
        raise NotImplementedError("This function has to be implemented in the derived class!")

    def ImportData(self, data_config):
        """Imports data from an external solver
        External solver sends, CoSimulation receives

        @param data_config <python dictionary> : configuration of the data to be imported
        """
        raise NotImplementedError("This function has to be implemented in the derived class!")

    def ExportData(self, data_config):
        """Exports data to an external solver
        CoSimulation sends, external solver receives

        @param data_config <python dictionary> : configuration of the data to be exported
        """
        raise NotImplementedError("This function has to be implemented in the derived class!")

    def PrintInfo(self):
        pass

    def Check(self):
        print("!!!WARNING!!! your IO does not implement Check!!!")

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"        : "",
            "echo_level"  : 0
        }""")
