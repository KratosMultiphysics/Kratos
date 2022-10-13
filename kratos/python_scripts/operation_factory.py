# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from importlib import import_module

class KratosOperationFactory(object):
    def __init__(self, model):
        self.model = model

    def ConstructListOfOperations( self, operations_list):
        constructed_operations = []

        for i in range(operations_list.size()):
            item = operations_list[i]

            if item.Has("operation"):
                """Location of the operations in Kratos; e.g.:
                - Registered c++ operation:                     operations.KratosMultiphysics.Operation
                - Registered c++ application operation:         operations.KratosMultiphysics.FluidDynamicsApplication.FluidOperation
                - Python operation:                             KratosMultiphysics.sample_operation.SampleOperation
                - Python application operation:                 KratosMultiphysics.FluidDynamicsApplication.fluid_sample_operation.FluidSampleOperation
                - Full path to custom operation:                user_operation.UserOperation (this has to be located in the simulation directory)
                """

                # Check if current operation is in the c++ registry
                # If not found, search for the corresponding Python module
                operation_string = item["operation"].GetString()
                operation_settings = item["Parameters"]
                if KM.Registry.HasItem(operation_string):
                    # operation_prototype = KM.Registry.GetOperation(operation_string) # TODO: Remove this after developing --> We want to access c++ registry through Python one
                    operation_prototype = KM.KratosGlobals.Registry.GetItem(operation_string)
                    operation = operation_prototype.Create(self.model, operation_settings)
                else:
                    operation_module_name, operation_name = operation_string.rsplit(".", 1)
                    try:
                        operation_module = import_module(operation_module_name)
                        operation = getattr(operation_module, operation_name).Create(self.model, operation_settings)
                    except ModuleNotFoundError as e:
                        raise ModuleNotFoundError(e)

                # Append the operation to the current operation list
                constructed_operations.append(operation)
            else:
                error_msg = f"Provided operation is incomplete: {item}.\n"
                error_msg += "'operation' field must be defined for each operation."
                raise NameError(error_msg)

        return constructed_operations
