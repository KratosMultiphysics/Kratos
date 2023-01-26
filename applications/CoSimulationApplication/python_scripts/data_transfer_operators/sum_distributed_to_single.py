# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Other imports
import numpy as np

def Create(*args):
    return SumDistributedToSingle(*args)

class SumDistributedToSingle(CoSimulationDataTransferOperator):
    """DataTransferOperator to sum values on one interface and put it to another interface.
    Used e.g. for FSI with SDof, where the loads on the fluid interface are summed up and set to the SDof interface
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        if not from_solver_data.IsDefinedOnThisRank():
            return

        data_array = from_solver_data.GetData()

        value = data_array.sum()
        value = from_solver_data.GetModelPart().GetCommunicator().GetDataCommunicator().Sum(value, 0)
        summed_data_array = np.array([value])

        if not to_solver_data.IsDefinedOnThisRank():
            return

        if to_solver_data.IsDistributed():
            raise Exception("The destination if the data cannot be distributed!")

        to_solver_data_size = to_solver_data.Size()
        if not to_solver_data_size == 1:
            raise Exception('Interface data "{}" of solver "{}" requires to be of size 1, got: {}'.format(to_solver_data.name, to_solver_data.solver_name, to_solver_data_size))

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            summed_data_array *= (-1)
        if "add_values" in transfer_options.GetStringArray():
            summed_data_array += to_solver_data.GetData()

        to_solver_data.SetData(summed_data_array)

    def _Check(self, from_solver_data, to_solver_data):
        if not from_solver_data.is_scalar_variable:
            raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(from_solver_data.name, from_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "add_values"]
