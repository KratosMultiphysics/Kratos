# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(*args):
    return CopyWithEmptyRanksDataTransferOperator(*args)

class CopyWithEmptyRanksDataTransferOperator(CoSimulationDataTransferOperator):
    """
    DataTransferOperator that copies values from one interface to another, without any checks.
    It is an extended version of CopyDataTransferOperation which considers solvers that only
    run in rank 0 (e.g. the RigidBodySolver).
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):

        # First, check if the sending solver has empty ranks
        # In this case, the data will be broadcasted from rankto all the 
        from_empty_rank = 0
        # If the rank is not empty, save the data as usual
        if from_solver_data.IsDefinedOnThisRank():
            from_solver_data_size = from_solver_data.Size()
            from_solver_data_array = from_solver_data.GetData()
        # If it is empty, put some dummy values to overwritte later
        else:
            from_empty_rank = 1
            from_solver_data_size = 0
            from_solver_data_array = []
        # Check if some ranks were empty. If not, the broadcast is not necessary
        empty_from_solver_ranks = self.data_communicator.SumAll(from_empty_rank)
        if empty_from_solver_ranks > 0:
            # Broadcast the data and its size
            from_solver_data_size = self.data_communicator.Broadcast(from_solver_data_size, 0)
            from_solver_data_array = self.data_communicator.BroadcastDoubles(from_solver_data_array, 0)

        # Second, check if the receiver solver has empty ranks so nothing is sent to them
        if not to_solver_data.IsDefinedOnThisRank():
            return

        # From here, it does the same as CopyDataTransferOperation
        to_solver_data_size = to_solver_data.Size()
        if not from_solver_data_size == to_solver_data_size:
            raise Exception('The sizes of the data are not matching: {} != {} for interface data "{}" of solver "{}" and interface data "{}" of solver "{}"!'.format(from_solver_data_size, to_solver_data_size, from_solver_data.name, from_solver_data.solver_name, to_solver_data.name, to_solver_data.solver_name))

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