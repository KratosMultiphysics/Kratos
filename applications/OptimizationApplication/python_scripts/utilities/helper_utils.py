from importlib import import_module
from pathlib import Path
from enum import Enum
from typing import Union

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class ContainerEnum(Enum):
    NODES = 1
    ELEMENTS = 2
    CONDITIONS = 3
    ELEMENT_PROPERTIES = 4
    CONDITION_PROPERTIES = 5

def GetSensitivityContainer(model_part: Kratos.ModelPart, container_type: ContainerEnum) -> 'Union[Kratos.ModelPart.Nodes, Kratos.ModelPart.Conditions, Kratos.ModelPart.Elements]':
    if container_type == ContainerEnum.NODES:
        return model_part.Nodes
    elif container_type in [ContainerEnum.ELEMENTS, ContainerEnum.ELEMENT_PROPERTIES]:
        return model_part.Elements
    elif container_type in [ContainerEnum.CONDITIONS, ContainerEnum.CONDITION_PROPERTIES]:
        return model_part.Conditions
    else:
        raise RuntimeError("Unsupported container type requested.")

class ContainerVariableDataHolder:
    def __init__(self, values: Kratos.Vector, variable, model_part: Kratos.ModelPart, container_type: ContainerEnum):
        self.__values = values
        self.__model_part = model_part
        self.__container_type = container_type
        self.__variable = variable

    def GetValues(self) -> Kratos.Vector:
        return self.__values

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.__model_part

    def GetContainerType(self) -> ContainerEnum:
        return self.__container_type

    def GetVariable(self):
        return self.__variable

    def GetContainer(self) -> 'Union[Kratos.ModelPart.Nodes, Kratos.ModelPart.Conditions, Kratos.ModelPart.Elements]':
        return GetSensitivityContainer(self.GetModelPart(), self.GetContainerType())

    def __str__(self) -> str:
        return f"ModelPart: {self.__model_part.FullName}, ContainerType: {self.__container_type.name}, VariableName: {self.__variable.Name()}"

def Factory(
    module_name: str,
    class_name: str,
    model: Kratos.Model,
    parameters: Kratos.Parameters,
    optimization_info: OptimizationInfo,
    required_object_type = object):

    python_file_name = ''.join(['_' + c.lower() if c.isupper() else c for c in class_name]).lstrip('_')
    full_module_name = f"{module_name}.{python_file_name}"

    # try importing the module
    try:
        module = import_module(full_module_name)
    except ModuleNotFoundError:
        RuntimeError(f"Invalid \"{module_name}\" and/or \"{class_name}\" provided with settings {str(parameters)} [ Full module name = {full_module_name} ].")

    # retrieving the object from the module
    try:
        retrieved_object = getattr(module, class_name)(model, parameters, optimization_info)
    except AttributeError:
        RuntimeError(f"The class \"{class_name}\" is not found in module \"{full_module_name}\"")

    module = import_module(full_module_name)
    retrieved_object = getattr(module, class_name)(model, parameters, optimization_info)

    # check retrieved object is of the required type
    if not isinstance(retrieved_object, required_object_type):
        raise RuntimeError(f"The retrieved object is of type \"{retrieved_object.__class__.__name__}\" which is not of the required type \"{required_object_type.__name__}\".")
    else:
        return retrieved_object

    # except AttributeError:
    #     raise AttributeError(f"The python file \"{str(Path(module.__file__))}\" does not have a class named \"{class_name}\".")
    # except (ImportError, ModuleNotFoundError):
    #     if module_prefix != "":
    #         try:
    #             module = import_module(module_prefix)
    #         except (ImportError, ModuleNotFoundError):
    #             raise RuntimeError(f"The module \"{module_prefix}\" cannot be imported.")
    #         if hasattr(module, "__path__"):
    #             list_of_available_options = []
    #             for item in Path(module.__path__[0]).iterdir():
    #                 if item.is_file() and item.name.endswith(".py"):
    #                     try:
    #                         current_module = import_module(module_prefix + "." + item.name[:-3])
    #                         current_class_name = "".join([c.title() for c in item.name[:-3].split("_")])
    #                         if hasattr(current_module, current_class_name):
    #                             if issubclass(getattr(current_module, current_class_name), required_object_type):
    #                                 list_of_available_options.append(current_class_name)
    #                     except:
    #                         pass
    #             raise RuntimeError(f"Given type with name \"{class_name}\" is not found. Followings are available options:\n\t" + "\n\t".join(list_of_available_options))
    #         else:
    #             raise RuntimeError(f"Invalid module found under \"{module_prefix}\".")
    #     else:
    #         raise RuntimeError(f"Invalid \"module\" and/or \"type\" fields. [ Provided settings = {parameters} ]")