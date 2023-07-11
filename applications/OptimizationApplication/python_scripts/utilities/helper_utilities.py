from pathlib import Path
from importlib import import_module
from typing import Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.kratos_utilities import GetListOfAvailableApplications
from KratosMultiphysics.kratos_utilities import GetKratosMultiphysicsPath
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

def GetClassModuleFromKratos(full_class_name: str) -> str:
    sub_module_paths = full_class_name.split(".")

    if not sub_module_paths:
        raise RuntimeError("Empty class names are not allowed.")

    relative_sub_module = Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(sub_module_paths[-1])
    if len(sub_module_paths) > 1:
        relative_sub_module = ".".join(sub_module_paths[:-1]) + f".{relative_sub_module}"

    relative_sub_module_path = relative_sub_module.replace(".", "/")
    kratos_path = GetKratosMultiphysicsPath()

    # check whether it is in Kratos core
    if Path(f"{kratos_path}/{relative_sub_module_path}.py").is_file():
        return f"KratosMultiphysics.{relative_sub_module}", sub_module_paths[-1]

    # now check if it is found in any of the compiled applications
    list_of_available_appliacations = GetListOfAvailableApplications()

    module_application = ""
    for application in list_of_available_appliacations:
        if Path(f"{kratos_path}/{application}/{relative_sub_module_path}.py").is_file():
            module_application = f"KratosMultiphysics.{application}.{relative_sub_module}"

    if module_application != "":
        return module_application, sub_module_paths[-1]
    else:
        raise RuntimeError(f"{full_class_name} is not found in KratosMultiphysics core or any of the available application directories in Kratos. Available applications:\n\t" + "\n\t".join(list_of_available_appliacations))

def CallOnAll(list_of_objects: 'list[Any]', method: Any, *args, **kwargs):
    for obj in list_of_objects:
        getattr(obj, method.__name__)(*args, **kwargs)

def IsSameContainerExpression(container_expression_1: ContainerExpressionTypes, container_expression_2: ContainerExpressionTypes) -> bool:
    if container_expression_1.GetModelPart().FullName() != container_expression_2.GetModelPart().FullName():
        return False

    if type(container_expression_1) != type(container_expression_2):
        return False

    return True

def HasContainerExpression(container_expression: ContainerExpressionTypes, list_of_container_expressions: 'list[ContainerExpressionTypes]') -> bool:
    return any([IsSameContainerExpression(container_expression, list_container_expression) for list_container_expression in list_of_container_expressions])

def OptimizationComponentFactory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    if not parameters.Has("type"):
        raise RuntimeError(f"Components created from OptimizationComponentFactory require the \"type\" [ provided paramters = {parameters}].")

    python_type = parameters["type"].GetString()

    if not parameters.Has("module") or parameters["module"].GetString() == "":
        # in the case python type comes without a module
        # as in the case python_type is in the sys path or the current working directory.
        full_module_name = python_type
    else:
        # in the case python type comes witha a module.
        module = parameters["module"].GetString()
        full_module_name = f"{module}.{python_type}"

    module = import_module(full_module_name)
    if not hasattr(module, "Factory"):
        raise RuntimeError(f"Python module {full_module_name} does not have a Factory method.")

    return getattr(module, "Factory")(model, parameters, optimization_problem)

def GetAllComponentFullNamesWithData(optimization_problem: OptimizationProblem) -> 'list[str]':
    data_container = optimization_problem.GetProblemDataContainer()

    list_of_components_full_names_with_data: 'list[str]' = []
    for component_type_str, components_dict in data_container.GetSubItems().items():
        if component_type_str != "object":
            component_type_str = Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(component_type_str) + "."
        else:
            component_type_str = ""
        for component_name in components_dict.GetSubItems().keys():
            list_of_components_full_names_with_data.append(f"{component_type_str}{component_name}")

    return list_of_components_full_names_with_data

def GetComponentHavingDataByFullName(component_full_name: str, optimization_problem: OptimizationProblem) -> Any:
    data_container = optimization_problem.GetProblemDataContainer()

    name_data = component_full_name.split(".")

    if len(name_data) == 1:
        if data_container.HasValue("object") and data_container["object"].HasValue(component_full_name):
            return component_full_name
    else:
        component_type_str = Kratos.StringUtilities.ConvertSnakeCaseToCamelCase(name_data[0])
        component_name = name_data[1]

        for component_type in optimization_problem.GetComponentContainer().keys():
            if component_type.__name__ == component_type_str:
                return optimization_problem.GetComponent(component_name, component_type)

    msg = ""
    for component_type, dict_of_components in optimization_problem.GetComponentContainer().items():
        for sub_item_name in dict_of_components.keys():
            msg += "\n\t" + Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(component_type.__name__) + f".{sub_item_name}"
    if data_container.HasValue("object"):
        for sub_item_name in data_container["object"].GetSubItems().keys():
            msg += "\n\t" + sub_item_name

    raise RuntimeError(f"\"{component_full_name}\" full component name is not found in the optimization problem. Followings are supported component with full names:" + msg)


