# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from importlib import import_module

################# Please do not change the following class:
class KratosProcessFactory(object):
    def __init__(self, Model):
        self.Model = Model #Model is a place

    def ConstructListOfProcesses( self, process_list ):
        constructed_processes = []
        for i in range(0,process_list.size()):
            item = process_list[i]

            # Registry-based instantiation
            if item.Has("name"):
                registry_entry = item["name"].GetString()
                if KM.Registry.HasItem(registry_entry):
                    # Get already stored prototype
                    if KM.Registry.HasItem(f"{registry_entry}.Prototype"):
                        prototype = KM.Registry[f"{registry_entry}.Prototype"]
                        instance = prototype.Create(self.Model, item["Parameters"])
                    # Get prototype from stored Python module
                    elif KM.Registry.HasItem(f"{registry_entry}.ModuleName"):
                        class_name = registry_entry.split(".")[-1]
                        module_name = KM.Registry[f"{registry_entry}.ModuleName"]
                        module = import_module(module_name)
                        if hasattr(module, class_name):
                            prototype = getattr(module, class_name)
                            instance = prototype(self.Model, item["Parameters"])
                        else:
                            #TODO: In here we're assuming that the registry last key is the class name
                            #TODO: We should enforce this. Now an error happens but as we populate we should throw a warning and search for a ClassName item
                            err_msg = f"The '{class_name}' class name cannot be found within the '{module_name}' module."
                            raise Exception(err_msg)
                    else:
                        err_msg = f"Registry process '{registry_entry}' cannot be constructed."
                        raise Exception(err_msg)
                    # Construct the process from the obtained prototype and append it to the list
                    constructed_processes.append(instance)
                else:
                    KM.Logger.PrintWarning(f"Asking to construct the non-registered 'name' '{registry_entry}'.")
            else:

                # Alternative (old) instantiation
                if not item.Has("python_module"):
                    KM.Logger.PrintWarning("Your list of processes: ", process_list)
                    raise NameError('"python_module" must be defined in your parameters. Check all your processes')

                # python-script that contains the process
                python_module_name = item["python_module"].GetString()

                if item.Has("kratos_module"): # for Kratos-processes
                    """Location of the process in Kratos; e.g.:
                    - KratosMultiphysics
                    - KratosMultiphysics.FluidDynamicsApplication
                    """

                    kratos_module_name = item["kratos_module"].GetString()
                    if not kratos_module_name.startswith("KratosMultiphysics"):
                        kratos_module_name = "KratosMultiphysics." + kratos_module_name

                    full_module_name = kratos_module_name + "." + python_module_name
                    python_module = import_module(full_module_name)

                    p = python_module.Factory(item, self.Model)
                    constructed_processes.append( p )

                else: # for user-defined processes
                    python_module = import_module(python_module_name)
                    p = python_module.Factory(item, self.Model)
                    constructed_processes.append( p )

        return constructed_processes


########## here we generate the common kratos processes --- IMPLEMENTED IN C++ ###################
def Factory(settings, Model):
    if(settings["process_name"].GetString() == "ApplyConstantScalarValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return KM.ApplyConstantScalarValueProcess(model_part, settings["Parameters"])

    elif(settings["process_name"].GetString() == "ApplyConstantVectorValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return KM.ApplyConstantVectorValueProcess(model_part, settings["Parameters"])
    elif(settings["process_name"].GetString() == "TimeAveragingProcess"):
        return KM.TimeAveragingProcess(Model, settings["Parameters"])

    raise Exception("Process name not found ",)
