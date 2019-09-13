from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

import numpy as np

def Create(settings):
    return SumDistToSingle(settings)

class SumDistToSingle(CoSimulationDataTransferOperator):

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        self._CheckAvailabilityTransferOptions(transfer_options)

        data_array = np.array([])
        data_array = from_solver_data.GetData()

        value = sum(data_array)
        if from_solver_data.IsDistributed():
            value = from_solver_data.GetModelPart.GetCommunicator().GetDataCommunicator().SumAll(value)
        summed_data_array = np.array([value])

        if "swap_sign" in transfer_options.GetStringArray():
            summed_data_array *= (-1)

        to_solver_data.SetData(summed_data_array)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign"]