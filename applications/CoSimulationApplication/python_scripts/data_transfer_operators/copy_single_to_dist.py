from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator
import numpy as np

def Create(settings):
    return CopySingleToDist(settings)

class CopySingleToDist(CoSimulationDataTransferOperator):

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        self._CheckAvailabilityTransferOptions(transfer_options)

        to_solver_values = to_solver_data.GetData()
        data_value = from_solver_data.GetData()

        if not data_value.size == 1:
            raise Exception('Expected one value, got: {}'.format(data_value.size))

        to_solver_values.fill(data_value[0])

        if "swap_sign" in transfer_options.GetStringArray():
            to_solver_values *= (-1)

        to_solver_data.SetData(to_solver_values)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign"]