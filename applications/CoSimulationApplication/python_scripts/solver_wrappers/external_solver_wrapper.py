from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_name):
    return ExternalSolverWrapper(settings, solver_name)

class ExternalSolverWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for external solvers
    """
    def __init__(self, settings, solver_name):
        super(ExternalSolverWrapper, self).__init__(settings, solver_name)

        settings_defaults = KM.Parameters("""{
            "model_parts_read" : { },
            "model_parts_send" : { },
            "model_parts_recv" : { }
        }""")

        self.settings["settings"].ValidateAndAssignDefaults(settings_defaults)

        self.__CreateModelPartsFromCouplingData()
        self._AllocateHistoricalVariablesFromCouplingData()

    def Initialize(self):
        for main_model_part_name, mdpa_file_name in self.settings["settings"]["model_parts_read"].items():
            KM.ModelPartIO(mdpa_file_name.GetString()).ReadModelPart(self.model[main_model_part_name])

        for model_part_name, comm_name in self.settings["settings"]["model_parts_send"].items():
            interface_config = {
                "comm_name" : comm_name.GetString(),
                "model_part_name" : model_part_name
            }

            self.ExportCouplingInterface(interface_config)

        for model_part_name, comm_name in self.settings["settings"]["model_parts_recv"].items():
            interface_config = {
                "comm_name" : comm_name.GetString(),
                "model_part_name" : model_part_name
            }

            self.ImportCouplingInterface(interface_config)

        super(ExternalSolverWrapper, self).Initialize()

    def AdvanceInTime(self, current_time):
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO

    def PrintInfo(self):
        cs_tools.cs_print_info(self._ClassName(), "printing info...")
        ## TODO print additional stuff with higher echo-level

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def __CreateModelPartsFromCouplingData(self):
        '''This function creates the Main-ModelParts that are used for external solvers
        '''
        for data in self.data_dict.values():
            main_model_part_name = data.model_part_name.split(".")[0]
            if not self.model.HasModelPart(main_model_part_name):
                self.model.CreateModelPart(main_model_part_name)
                if self.echo_level > 0:
                    cs_tools.cs_print_info("ExternalSolverWrapper", 'Created ModelPart "{}" for solver "{}"'.format(main_model_part_name, self.name))
