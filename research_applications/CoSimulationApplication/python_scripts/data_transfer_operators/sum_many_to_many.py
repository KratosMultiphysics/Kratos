# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(*args):
    return SumManyToMany(*args)

class SumManyToMany(CoSimulationDataTransferOperator):
    """DataTransferOperator to sum values on one interface and put it to another interface, if desired equally distributed.    
    """

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        if not from_solver_data.IsDefinedOnThisRank():
            return

        data_values = from_solver_data.GetData()
        data_value = data_values.sum()        

        if from_solver_data.IsDistributed():
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

        # check the to_solver_data
        if to_solver_data.IsDefinedOnThisRank():
            if not to_solver_data.is_scalar_variable:
                raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(to_solver_data.name, to_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "distribute_values", "add_values"]