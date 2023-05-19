from pathlib import Path
from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.kratos_utilities import GetListOfAvailableApplications
from KratosMultiphysics.kratos_utilities import GetKratosMultiphysicsPath
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

def GetClassModuleFromKratos(class_name: str) -> str:
    snake_case_class_name = Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(class_name)
    kratos_path = GetKratosMultiphysicsPath()

    # check whether it is in Kratos core
    if Path(f"{kratos_path}/{snake_case_class_name}.py").is_file():
        return "KratosMultiphysics"

    # now check if it is found in any of the compiled applications
    list_of_available_appliacations = GetListOfAvailableApplications()

    module_application = ""
    for application in list_of_available_appliacations:
        if Path(f"{kratos_path}/{application}/{snake_case_class_name}.py").is_file():
            module_application = f"KratosMultiphysics.{application}"

    if module_application != "":
        return module_application
    else:
        raise RuntimeError(f"{class_name} is not found in KratosMultiphysics core or any of the available application root directories in Kratos. Available applications:\n\t" + "\n\t".join(list_of_available_appliacations))

def CallOnAll(list_of_objects: 'list[any]', method: any, *args, **kwargs):
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