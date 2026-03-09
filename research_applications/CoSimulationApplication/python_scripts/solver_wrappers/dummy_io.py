# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

def Create(*args):
    return DummyIO(*args)

class DummyIO(CoSimulationIO):
    """This class is used if a Solver directly uses Kratos as a data-structure
    e.g. Kratos itself or simple-solvers written in Python
    """

    def ImportCouplingInterface(self, interface_config):
        pass

    def ExportCouplingInterface(self, interface_config):
        pass

    def ImportData(self, data_config):
        pass

    def ExportData(self, data_config):
        pass

    def PrintInfo(self):
        print("This is the dummy-IO")

    def Check(self):
        pass
