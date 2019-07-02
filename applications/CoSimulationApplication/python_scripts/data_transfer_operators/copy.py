from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(settings):
    return CopyDataTransferOperator(settings)

class CopyDataTransferOperator(CoSimulationDataTransferOperator):
    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        # TODO to be implemented, depending on the location of the data
        raise NotImplementedError
