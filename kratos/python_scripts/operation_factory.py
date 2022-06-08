# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from importlib import import_module

################# Please do not change the following class:
class KratosOperationFactory(object):
    def __init__(self, model):
        self.model = model #Model is a beautiful place

    def ConstructListOfOperations( self, operations_list):
        constructed_operations = []
        for i in range(operations_list.size()):
            item = operations_list[i]

            if not item.Has("python_module"):
                pass
            else:
                # python-script that contains the process
                python_module_name = item["python_module"].GetString()

                if item.Has("kratos_module"): # for Kratos-operations
                    """Location of the operations in Kratos; e.g.:
                    - KratosMultiphysics
                    - KratosMultiphysics.FluidDynamicsApplication
                    """

                    kratos_module_name = item["kratos_module"].GetString()
                    if not kratos_module_name.startswith("KratosMultiphysics"):
                        kratos_module_name = "KratosMultiphysics." + kratos_module_name

                    full_module_name = kratos_module_name + "." + python_module_name
                    python_module = import_module(full_module_name)

                    p = python_module.Factory(item, self.model)
                    constructed_operations.append( p )

                else: # for user-defined operations
                    python_module = import_module(python_module_name)
                    p = python_module.Factory(item, self.model)
                    constructed_operations.append( p )

        return constructed_operations
