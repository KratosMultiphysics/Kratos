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

        self.__CreateModelPartsFromCouplingData()
        self._AllocateHistoricalVariablesFromCouplingData()

    def Initialize(self):
        # receive coupling-interface from external solvers
        for data in self.data_dict.values():
            self.io.ImportCouplingInterface(data.GetModelPart())

        super(ExternalSolverWrapper, self).Initialize()

    def AdvanceInTime(self, current_time):
        return 0.0

    def SolveSolutionStep(self):
        # Send data-fields to external solver
        for data in self.data_dict.values():
            self.io.ExportCouplingInterfaceData(data)

        # External Solver solves

        # Receive data-fields from external solver
        for data in self.data_dict.values():
            self.io.ImportCouplingInterfaceData(data)

    def __CreateModelPartsFromCouplingData(self):
        '''This function creates the Main-ModelParts that are used for external solvers
        '''
        for data in self.data_dict.values():
            main_model_part_name = data.model_part_name.split(".")[0]
            if not self.model.HasModelPart(main_model_part_name):
                self.model.CreateModelPart(main_model_part_name)
                if self.echo_level > 0:
                    cs_tools.cs_print_info("ExternalSolverWrapper", 'Created ModelPart "{}" for solver "{}"'.format(main_model_part_name, self.name))
