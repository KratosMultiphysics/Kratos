# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

import pdb

def Create(model, settings, solver_name):
    return PfemIO(model, settings, solver_name)

class PfemIO(CoSimulationIO):
    """This class is used if a Solver directly uses Kratos as a data-structure
    e.g. Kratos itself or simple-solvers written in Python
    """

    def ImportCouplingInterface(self, interface_config):
        pass

    def ExportCouplingInterface(self, interface_config):
        pass

    def ImportData(self, data_config):
        pdb.set_trace()
        new_data_config = data_config
        new_data_config = self.repeat_time_step_flag
        return new_data_config

    def ExportData(self, data_config):
        #pdb.set_trace()
        self.repeat_time_step_flag = data_config["repeat_time_step"]

    def PrintInfo(self):
        print("This is the pfem-IO")

    def Check(self):
        pass
