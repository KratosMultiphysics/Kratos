# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

def Create(*args):
    return ExternalPythonSolverIO(*args)

class ExternalPythonSolverIO(CoSimulationIO):
    """This class is used if a Solver directly uses Kratos as a data-structure
    e.g. Kratos itself or simple-solvers written in Python
    """
    def __init__(self, settings, model, solver_name, data_communicator):
        super().__init__(settings, model, solver_name, data_communicator)
    def ImportCouplingInterface(self, interface_config):
        pass

    def ExportCouplingInterface(self, interface_config):
        pass

    def ImportData(self, data_config):
        external_model = data_config["external_model"]
        interface_data = data_config["interface_data"]
        external_model.ImportData(interface_data)
        # raise RuntimeError("ImportData")

    def ExportData(self, data_config):
        if data_config["type"] == "direct_value":
            external_model = data_config["external_model"]
            interface_data = data_config["interface_data"]
            external_model.ExportData(interface_data)

    def PrintInfo(self):
        print("This is the dummy-IO")

    def Check(self):
        pass
