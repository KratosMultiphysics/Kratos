from importlib import import_module
from pathlib import Path
from enum import Enum

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class ContainerEnum(Enum):
    NODES = 1
    ELEMENTS = 2
    CONDITIONS = 3
    ELEMENT_PROPERTIES = 4
    CONDITION_PROPERTIES = 5

def GetSensitivityContainer(model_part: Kratos.ModelPart, container_type: ContainerEnum):
    match (container_type):
        case ContainerEnum.NODES:
            return model_part.Nodes
        case ContainerEnum.ELEMENTS | ContainerEnum.ELEMENT_PROPERTIES:
            return model_part.Elements
        case ContainerEnum.CONDITIONS | ContainerEnum.CONDITION_PROPERTIES:
            return model_part.Conditions
        case _:
            raise RuntimeError("Unsupported container type requested.")

def RetrieveObject(model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo, object_type = object):
    default_settings = Kratos.Parameters("""{
        "module"  : "",
        "type"    : "",
        "settings": {}
    }""")
    parameters.ValidateAndAssignDefaults(default_settings)

    class_name = parameters["type"].GetString()
    python_file_name = ''.join(['_' + c.lower() if c.isupper() else c for c in class_name]).lstrip('_')

    module_prefix = parameters["module"].GetString()
    module_name = module_prefix
    if module_name != "":
        module_name += "."
    module_name += python_file_name

    try:
        module = import_module(module_name)
        retrieved_object = getattr(module, class_name)(model, parameters["settings"], optimization_info)
        if not isinstance(retrieved_object, object_type):
            raise ImportError(f"The retrieved object is of type \"{retrieved_object.__class__.__name__}\" which is not of the required type \"{object_type.__name__}\".")
        else:
            return retrieved_object
    except AttributeError:
        raise AttributeError(f"The python file \"{str(Path(module.__file__))}\" does not have a class named \"{class_name}\".")
    except (ImportError, ModuleNotFoundError):
        if module_prefix != "":
            try:
                module = import_module(module_prefix)
            except (ImportError, ModuleNotFoundError):
                raise RuntimeError(f"The module \"{module_prefix}\" cannot be imported.")
            if hasattr(module, "__path__"):
                list_of_available_options = []
                for item in Path(module.__path__[0]).iterdir():
                    if item.is_file() and item.name.endswith(".py"):
                        try:
                            current_module = import_module(module_prefix + "." + item.name[:-3])
                            current_class_name = "".join([c.title() for c in item.name[:-3].split("_")])
                            if hasattr(current_module, current_class_name):
                                if issubclass(getattr(current_module, current_class_name), object_type):
                                    list_of_available_options.append(current_class_name)
                        except:
                            pass
                raise RuntimeError(f"Given type with name \"{class_name}\" is not found. Followings are available options:\n\t" + "\n\t".join(list_of_available_options))
            else:
                raise RuntimeError(f"Invalid module found under \"{module_prefix}\".")
        else:
            raise RuntimeError("Invalid \"module\" and/or \"type\" fields.")