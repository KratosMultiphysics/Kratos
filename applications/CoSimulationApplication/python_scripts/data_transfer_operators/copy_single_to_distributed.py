# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(*args):
    return CopySingleToDistributed(*args)

class CopySingleToDistributed(CoSimulationDataTransferOperator):
    """DataTransferOperator to take one single value and set it to all values on another interface.
    Used e.g. for FSI with SDof, where the SDof has one value and the fluid interface has many
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        if not to_solver_data.IsDefinedOnThisRank():
            return

        if from_solver_data.IsDefinedOnThisRank():
            data_value = from_solver_data.GetData()

            if from_solver_data.IsDistributed():
                raise Exception("The source if the data cannot be distributed!")

            if not data_value.size == 1:
            # assert from_solver_data.Size() == 1
                raise Exception('Interface data "{}" of solver "{}" requires to be of size 1, got: {}'.format(from_solver_data.name, from_solver_data.solver_name, data_value.size))
            data_value = data_value[0]
        else:
            data_value = 0.0

        data_comm = to_solver_data.model_part.GetCommunicator().GetDataCommunicator()
        data_value = data_comm.Broadcast(data_value, 0)

        to_solver_values = to_solver_data.GetData()
        to_solver_values.fill(data_value)

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            to_solver_values *= (-1)
        if "distribute_values" in transfer_options.GetStringArray():
            to_solver_values /= data_comm.SumAll(to_solver_data.Size())
        if "add_values" in transfer_options.GetStringArray():
            to_solver_values += to_solver_data.GetData()

        to_solver_data.SetData(to_solver_values)

    def _Check(self, from_solver_data, to_solver_data):
        if not to_solver_data.IsDefinedOnThisRank():
            return
        if not to_solver_data.is_scalar_variable:
            raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(to_solver_data.name, to_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "distribute_values", "add_values"]
