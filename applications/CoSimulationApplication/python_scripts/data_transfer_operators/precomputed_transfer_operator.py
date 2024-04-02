import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

def Create(*args):
    return PrecomputedTransferOperator(*args)

class PrecomputedTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that uses precomputed transfer matrices to transfer data between interfaces.
    """
    def __init__(self, settings, parent_coupled_solver_data_communicator):
        super().__init__(settings, parent_coupled_solver_data_communicator)

        # Initialize a dictionary to store the transfer operators
        self.transfer_operators = {}

        # Load the transfer operator matrices from the specified files
        if settings.Has("transfer_operator_files"):
            transfer_operator_files = settings["transfer_operator_files"]
            for variable_name in transfer_operator_files.GetKeys():
                file_path = transfer_operator_files[variable_name].GetString()
                self.transfer_operators[variable_name] = np.load(file_path)
        else:
            raise Exception('No "transfer_operator_files" provided!')

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        variable_name = from_solver_data.variable.Name()

        # Check if the transfer operator for the given variable exists
        if variable_name not in self.transfer_operators:
            raise Exception(f"No transfer operator available for variable {variable_name}")

        # Retrieve the data as a NumPy array
        data = from_solver_data.GetValues()

        # Perform the transfer operation using the precomputed matrix
        transferred_data = np.dot(self.transfer_operators[variable_name], data)

        # Set the computed values in the destination
        to_solver_data.SetValues(transferred_data, transfer_options)

    def _Check(self, from_solver_data, to_solver_data):
        # Add any checks necessary for your specific use case, e.g., matching dimensions
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        # Define and return the default parameters for this operator, including multiple transfer operators
        this_defaults = KM.Parameters("""{
            "transfer_operator_files": {
                "TEMPERATURE": "path/to/temperature_transfer_operator.npy",
                "AUX_FLUX": "path/to/flux_transfer_operator.npy"
            }
        }""")
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return cls.__mapper_flags_dict.keys()

    def __GetMapperFlags(self, transfer_options, from_solver_data, to_solver_data):
        mapper_flags = KM.Flags()
        for flag_name in transfer_options.GetStringArray():
            mapper_flags |= self.__mapper_flags_dict[flag_name]
        if from_solver_data.location == "node_non_historical":
            mapper_flags |= KM.Mapper.FROM_NON_HISTORICAL
        if to_solver_data.location == "node_non_historical":
            mapper_flags |= KM.Mapper.TO_NON_HISTORICAL

        return mapper_flags

    @staticmethod
    def __GetModelPartFromInterfaceData(interface_data):
        # Directly return the ModelPart from the interface data
        return interface_data.GetModelPart()
