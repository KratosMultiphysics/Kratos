from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator
import numpy as np

def Create(settings):
    return CopySingleToDist(settings)

class CopySingleToDist(CoSimulationDataTransferOperator):

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        self._CheckAvailabilityTransferOptions(transfer_options)

        size = to_solver_data.Size()
        transfer_options_list = transfer_options.GetStringArray()
        data_value = from_solver_data.GetData()

        data_array = np.zeros(size)
        for i in range(0, size, 1):
            data_array[i] = data_value[0]

        if "swap_sign" in transfer_options_list:
            data_value *= (-1)
            to_solver_data.SetData(data_array)
        else:
            to_solver_data.SetData(data_array)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign"]