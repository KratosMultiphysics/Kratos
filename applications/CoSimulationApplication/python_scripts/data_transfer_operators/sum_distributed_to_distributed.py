# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Other imports
import numpy as np

def Create(settings):
    return SumDistributedToDistributed(settings)

class SumDistributedToDistributed(CoSimulationDataTransferOperator):
    """DataTransferOperator to sum values on one interface and put equally distributed to another interface.    
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):        

        data_array = from_solver_data.GetData()
        to_solver_values = to_solver_data.GetData()

        #Sum up values from model_part from all ranks
        value = float(sum(data_array))
        if from_solver_data.IsDistributed():
            value = from_solver_data.model_part.GetCommunicator().GetDataCommunicator().SumAll(value)

        # distribute to new model_part
        to_solver_values.fill(value)

        # distribute equally on all nodes from model_part from all ranks
        data_comm = from_solver_data.model_part.GetCommunicator().GetDataCommunicator()
        to_solver_values /= data_comm.SumAll(to_solver_data.Size())

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            to_solver_values *= (-1)
        if "add_values" in transfer_options.GetStringArray():
            to_solver_values += to_solver_data.GetData()

        to_solver_data.SetData(to_solver_values)

    def _Check(self, from_solver_data, to_solver_data):
        if not from_solver_data.is_scalar_variable:
            raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(from_solver_data.name, from_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "add_values"]
