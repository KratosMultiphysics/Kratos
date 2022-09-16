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
                - Full path in Kratos:              KratosMultiphysics.operations.sample_operation.SampleOperation
                - Full path to custom operation:    user_operation.UserOperation  (this has to be located in your simulation directory)
                """

                operation_string = item["operation"].GetString()
                operation_module, operation_name = operation_string.rsplit(".", 1)
                print("Importing operation:", operation_module, operation_name)

                module = import_module(operation_module)

                p = module.Create(operation_name, self.model)
                constructed_operations.append( p )
            else:
                KM.Logger.PrintWarning("Your list of operations: ", operations_list)
                raise NameError('"operation" must be defined in your parameters. Check all your operations')

        return constructed_operations
