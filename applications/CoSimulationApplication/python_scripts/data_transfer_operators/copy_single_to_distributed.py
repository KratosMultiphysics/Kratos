from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

def Create(settings):
    return CopySingleToDistributed(settings)

class CopySingleToDistributed(CoSimulationDataTransferOperator):
    """DataTransferOperator to take one single value and set it to all values on another interface.
    Used e.g. for FSI with SDof, where the SDof has one value and the fluid interface has many
    """
    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        to_solver_values = to_solver_data.GetData()
        data_value = from_solver_data.GetData()

        if not data_value.size == 1:
            raise Exception('Interface data "{}" of solver "{}" requires to be of size 1, got: {}'.format(from_solver_data.name, from_solver_data.solver_name, data_value.size))

        to_solver_values.fill(data_value[0])

        # the order is IMPORTANT here!
        if "swap_sign" in transfer_options.GetStringArray():
            to_solver_values *= (-1)
        if "distribute_values" in transfer_options.GetStringArray():
            data_comm = from_solver_data.model_part.GetCommunicator().GetDataCommunicator()
            to_solver_values /= data_comm.SumAll(to_solver_data.Size())
        if "add_values" in transfer_options.GetStringArray():
            to_solver_values += to_solver_data.GetData()

        to_solver_data.SetData(to_solver_values)

    def _Check(self, from_solver_data, to_solver_data):
        if not to_solver_data.is_scalar_variable:
            raise Exception('Variable of interface data "{}" of solver "{}" has to be a scalar!'.format(to_solver_data.name, to_solver_data.solver_name))

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign", "distribute_values", "add_values"]
