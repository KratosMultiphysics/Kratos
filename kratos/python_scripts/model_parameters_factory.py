# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from importlib import import_module

class KratosModelParametersFactory(object):
    ''' Registry-based model and Kratos parameters factory
    This class implements a registry-based factory for any model
    and parameters constructible object (e.g. processes, operations)
    '''

    def __init__(self, model):
        self.model = model

    def ConstructItem(self, item):
        # Check the input settings
        # Note that we support both 'parameters' and 'Parameters' to keep processes backward compatibility
        if item.Has("parameters"):
            parameters = item["parameters"]
        elif item.Has("Parameters"):
            #TODO: Eventually throw a deprecation warning in here
            parameters = item["Parameters"]
        else:
            error_msg = f"Provided item '{item}' is incomplete.\n"
            error_msg += "'parameters' field must be defined for each item."
            raise NameError(error_msg)

        # Check for the 'name' entry
        # Note that this is a registry-based factory so 'name' must be always provided
        if item.Has("name"):
            """Location of the item in Kratos; e.g.:
            - Registered operation:                     operations.KratosMultiphysics.Operation
            - Registered application process:           processes.KratosMultiphysics.FluidDynamicsApplication.FluidProcess
            - Full path to Kratos modeler:              KratosMultiphysics.modelers.delete_model_parts_modeler.DeleteModelPartsModeler
            - Full path to custom operation:            user_operation.UserOperation (this has to be located in the simulation directory)
            """

            # Check if current item is registered
            # If not found, search for the corresponding Python module
            registry_entry = item["name"].GetString()
            if KM.Registry.HasItem(registry_entry):
                # Get already stored prototype
                if KM.Registry.HasItem(f"{registry_entry}.Prototype"):
                    prototype = KM.Registry[f"{registry_entry}.Prototype"]
                    instance = prototype.Create(self.model, parameters)
                # Get the create function
                elif KM.Registry.HasItem(f"{registry_entry}.CreateFunction"):
                    create = KM.Registry[f"{registry_entry}.CreateFunction"]
                    instance = create(self.model, parameters)
                # Get prototype from stored Python module
                elif KM.Registry.HasItem(f"{registry_entry}.ModuleName"):
                    class_name = registry_entry.split(".")[-1]
                    module_name = KM.Registry[f"{registry_entry}.ModuleName"]
                    module = import_module(module_name)
                    # Check if the standard class name is available
                    # Note that in here we assume that the registry last key matches the class name
                    if hasattr(module, class_name):
                        instance = getattr(module, class_name)(self.model, parameters)
                    else:
                        # If the registry key does not contain the class name we check for the 'ClassName' entry
                        if KM.Registry.HasItem(f"{registry_entry}.ClassName"):
                            class_name = KM.Registry[f"{registry_entry}.ClassName"]
                            KM.Logger.PrintWarning(f"Trying to construct item from non-standard 'ClassName' value {class_name}.")
                            instance = getattr(module, class_name)(self.model, parameters)
                        else:
                            err_msg = f"The '{class_name}' class name cannot be found within the '{module_name}' module."
                            raise Exception(err_msg)
            else:
                item_module_name, item_class_name = registry_entry.rsplit(".", 1)
                try:
                    item_module = import_module(item_module_name)
                    instance = getattr(item_module, item_class_name)(self.model, parameters)
                except ModuleNotFoundError as e:
                    raise ModuleNotFoundError(e)

            return instance
        else:
            error_msg = f"Provided item '{item}' is incomplete.\n"
            error_msg += "'name' field must be defined for each item."
            raise NameError(error_msg)

    def ConstructListOfItems( self, items_list):
        constructed_items = []

        for item in items_list:
            constructed_items.append(self.ConstructItem(item))

        return constructed_items
