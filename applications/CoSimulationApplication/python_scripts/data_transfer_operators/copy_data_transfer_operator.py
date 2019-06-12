from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_data_transfer_operator import CoSimulationBaseDataTransferOperator

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(settings):
    return CopyDataTransferOperator(settings)

class CopyDataTransferOperator(CoSimulationBaseDataTransferOperator):
    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        # TODO to be implemented, depending on the location of the data
        pass
        raise NotImplementedError
