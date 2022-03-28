# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(*args):
    return TransferOneToMany(*args)

class TransferOneToMany(CoSimulationDataTransferOperator):
    """DataTransferOperator to take one single value and set it to all values on another interface.
    Used e.g. for FSI with SDof, where the SDof has one value and the fluid interface has many
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        data_value = 0.0
        if from_solver_data.IsDefinedOnThisRank():
            data_values = from_solver_data.GetData()
            if data_values.size == 1: # this is the rank that actually contains the value
                data_value = data_values[0]

        data_value = self.data_communicator.SumAll(data_value)

        if not to_solver_data.IsDefinedOnThisRank():
            return

        to_solver_values = to_solver_data.GetData()
        to_solver_values.fill(data_value)

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            to_solver_values *= (-1)
        if "distribute_values" in transfer_options.GetStringArray():
            to_solver_values /= self.data_communicator.SumAll(to_solver_data.Size())
        if "add_values" in transfer_options.GetStringArray():
            to_solver_values += to_solver_data.GetData()

        to_solver_data.SetData(to_solver_values)

    def _Check(self, from_solver_data, to_solver_data):
        # check the from_solver_data
        if from_solver_data.IsDefinedOnThisRank():
            if not from_solver_data.is_scalar_variable:
                raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(from_solver_data.name, from_solver_data.solver_name))

        # make sure there is only one value
        if from_solver_data.IsDefinedOnThisRank():
            data_size = from_solver_data.Size()
            if from_solver_data.IsDistributed():
                data_size = from_solver_data.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(data_size)

            if not data_size == 1:
                raise Exception('Interface data "{}" of solver "{}" requires to be of size 1, got: {}'.format(from_solver_data.name, from_solver_data.solver_name, data_size))

        # check the to_solver_data
        if to_solver_data.IsDefinedOnThisRank():
            if not to_solver_data.is_scalar_variable:
                raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(to_solver_data.name, to_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "distribute_values", "add_values"]
