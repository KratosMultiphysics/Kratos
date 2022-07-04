# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Other imports
import numpy as np

def Create(settings):
    return SumDistributedToSingleVectorial(settings)

class SumDistributedToSingleVectorial(CoSimulationDataTransferOperator):
    """DataTransferOperator to sum values on one interface and put it to another interface.
    It is a version of SumDistributedToSingle which does the same but for vectorial variables.
    Used e.g. for FSI with RigidBody, where the loads on the fluid interface are summed up
    and set to the RigidBody interface
    """
    
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):

        data_array = from_solver_data.GetData().reshape(-1,from_solver_data.dimension)

        summed_data_array = []
        for dim in range(from_solver_data.dimension):

            value = sum(data_array[:,dim])
            if from_solver_data.IsDistributed():
                value = from_solver_data.GetModelPart.GetCommunicator().GetDataCommunicator().SumAll(value)
            dim_summed_data_array = np.array([value])

            # the order is IMPORTANT here!
            if "swap_sign" in transfer_options.GetStringArray():
                dim_summed_data_array *= (-1)
            if "add_values" in transfer_options.GetStringArray():
                dim_summed_data_array += to_solver_data.GetData()

            summed_data_array.append(dim_summed_data_array)

        to_solver_data.SetData(summed_data_array)

    def _Check(self, from_solver_data, to_solver_data):
        if to_solver_data.dimension != from_solver_data.dimension:
            msg = 'Interface data "' + str(to_solver_data.name) + '" of solver "' + str(to_solver_data.solver_name) + '" requires to be '
            msg += 'the same size as interface data "' + str(from_solver_data.name) + '" of solver "' + str(from_solver_data.solver_name)
            msg += '", got ' + str(to_solver_data.dimension) + ' and ' + str(from_solver_data.dimension) + ' respectively.'
            raise Exception(msg)

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "add_values"]