from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from ..base_classes.co_simulation_io import CoSimulationIO

def Create(model, custom_settings):
    return DummyIO(model, custom_settings)

class DummyIO(CoSimulationIO):
    """This class is used if a Solver directly uses Kratos as a data-structure
    e.g. Kratos itself or simple-solvers written in Python
    """

    def ImportCouplingInterfaceData(self, data_object, from_solver=None):
        pass

    def ImportCouplingInterface(self, mesh_config, from_solver=None):
        pass

    def ExportCouplingInterfaceData(self, data_object, to_solver=None):
        pass

    def ExportCouplingInterface(self, mesh_config, to_solver=None):
        pass

    def PrintInfo(self):
        print("This is the dummy-IO")

    def Check(self):
        pass
