# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(settings):
    return CopySingleToDistributedVectorial(settings)

class CopySingleToDistributedVectorial(CoSimulationDataTransferOperator):
    """DataTransferOperator to take one single value and set it to all values on another interface.
    It is a version of CopySingleToDistributed which does the same but for vectorial variables.
    Used e.g. for FSI with RigidBody, where the RigidBody has one value and the fluid interface has many.
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        to_solver_values = to_solver_data.GetData()
        data_value = from_solver_data.GetData()

        data_comm = from_solver_data.model_part.GetCommunicator().GetDataCommunicator()

        if data_value.size == 0:
            value = [float(0)] * from_solver_data.dimension
        else:
            value = data_value

        # sum Up Value from solver from all ranks -> so not zero when not definied in this rank
        value = from_solver_data.model_part.GetCommunicator().GetDataCommunicator().SumAll(value)

        for to_index in range(len(to_solver_values)):
            from_index = to_index % from_solver_data.dimension
            to_solver_values[to_index] = value[from_index]

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            to_solver_values *= (-1)
        if "distribute_values" in transfer_options.GetStringArray():
            to_solver_values /= (data_comm.SumAll(to_solver_data.Size()) / to_solver_data.dimension)
        if "add_values" in transfer_options.GetStringArray():
            to_solver_values += to_solver_data.GetData()

        to_solver_data.SetData(to_solver_values)

    def _Check(self, from_solver_data, to_solver_data):
        if to_solver_data.dimension != from_solver_data.dimension:
            msg = 'Interface data "' + str(to_solver_data.name) + '" of solver "' + str(to_solver_data.solver_name) + '" requires to be '
            msg += 'the same size as interface data "' + str(from_solver_data.name) + '" of solver "' + str(from_solver_data.solver_name)
            msg += '", got ' + str(to_solver_data.dimension) + ' and ' + str(from_solver_data.dimension) + ' respectively.'
            raise Exception(msg)
        if from_solver_data.Size() != from_solver_data.dimension:
            msg = 'Variable of interface data "{}" of solver "{}" has to be as large as its dimension.'
            raise Exception(msg.format(from_solver_data.name, from_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "distribute_values", "add_values"]
