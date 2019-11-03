from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Other imports
import numpy as np

def Create(settings):
    return SumDistributedToSingle(settings)

class SumDistributedToSingle(CoSimulationDataTransferOperator):
    """DataTransferOperator to sum values on one interface and put it to another interface.
    Used e.g. for FSI with SDof, where the loads on teh fluid interface are summed up and set to the SDofinterface has many
    """
    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        self._CheckAvailabilityTransferOptions(transfer_options)

        data_array = np.array([])
        data_array = from_solver_data.GetData()

        value = sum(data_array)
        if from_solver_data.IsDistributed():
            value = from_solver_data.GetModelPart.GetCommunicator().GetDataCommunicator().SumAll(value)
        summed_data_array = np.array([value])

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            summed_data_array *= (-1)
        if "add_values" in transfer_options.GetStringArray():
            summed_data_array += to_solver_data.GetData()

        to_solver_data.SetData(summed_data_array)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "add_values"]