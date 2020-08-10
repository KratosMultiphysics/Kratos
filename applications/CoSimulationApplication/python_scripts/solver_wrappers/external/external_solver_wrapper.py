from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, model, solver_name):
    return ExternalSolverWrapper(settings, model, solver_name)

class ExternalSolverWrapper(CoSimulationSolverWrapper):
    """This class is a generic wrapper for connecting external solvers
    The import of meshes is done once in the beginning
    """
    def __init__(self, settings, model, solver_name):
        super(ExternalSolverWrapper, self).__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "import_meshes"    : [ ]
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)

        cs_tools.CreateModelPartsFromCouplingData(self.data_dict.values(), self.model, self.name)
        cs_tools.AllocateHistoricalVariablesFromCouplingData(self.data_dict.values(), self.model, self.name)

    def Initialize(self):
        super(ExternalSolverWrapper, self).Initialize()

        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = { "model_part_name" : model_part_name }
            self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO

    def PrintInfo(self):
        cs_tools.cs_print_info(self._ClassName(), "printing info...")
        ## TODO print additional stuff with higher echo-level

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()