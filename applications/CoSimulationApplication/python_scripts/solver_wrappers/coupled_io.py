# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

def Create(model, settings, solver_name):
    return CoupledIO(model, settings, solver_name)

class CoupledIO(CoSimulationIO):
    """This class is used as an IO for any of the coupled solvers. The tricky part of accessing
        the data from involved solvers implemented here
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
        print("This is the Coupled-IO")

    def Check(self):
        pass
