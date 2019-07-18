from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(settings):
    return CopyDataTransferOperator(settings)

class CopyDataTransferOperator(CoSimulationDataTransferOperator):
    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        self._CheckAvailabilityTransferOptions(transfer_options)

        from_solver_data_size = from_solver_data.Size()
        to_solver_data_size = to_solver_data.Size()
        if not from_solver_data_size == to_solver_data_size:
            raise Exception('The sizes of the data are not matching: {} != {}!'.format(from_solver_data_size, to_solver_data_size))

        from_solver_data_array = from_solver_data.GetData()

        transfer_options_list = transfer_options.GetStringArray()

        if "swap_sign" in transfer_options_list:
            from_solver_data_array *= (-1)

        if "add_values" in transfer_options_list:
            to_solver_data.SetData(to_solver_data.GetData() + from_solver_data_array)
        else:
            to_solver_data.SetData(from_solver_data_array)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "add_values"]