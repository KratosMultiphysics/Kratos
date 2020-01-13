from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(settings):
    return CopyDataTransferOperator(settings)

class CopyDataTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that copies values from one interface to another, without any checks
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        from_solver_data_size = from_solver_data.Size()
        to_solver_data_size = to_solver_data.Size()
        if not from_solver_data_size == to_solver_data_size:
            raise Exception('The sizes of the data are not matching: {} != {} for interface data "{}" of solver "{}" and interface data "{}" of solver "{}"!'.format(from_solver_data_size, to_solver_data_size, from_solver_data.name, from_solver_data.solver_name, to_solver_data.name, to_solver_data.solver_name))

        from_solver_data_array = from_solver_data.GetData()

        transfer_options_list = transfer_options.GetStringArray()

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options_list:
            from_solver_data_array *= (-1)
        if "add_values" in transfer_options.GetStringArray():
            from_solver_data_array += to_solver_data.GetData()

        to_solver_data.SetData(from_solver_data_array)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "add_values"]