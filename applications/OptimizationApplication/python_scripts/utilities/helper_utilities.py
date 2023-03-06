from pathlib import Path
from importlib import import_module

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.kratos_utilities import GetListOfAvailableApplications
from KratosMultiphysics.kratos_utilities import GetKratosMultiphysicsPath

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

def OptimizationProcessFactory(
    module_name: str,
    class_name: str,
    model: Kratos.Model,
    parameters: Kratos.Parameters,
    optimization_info: OptimizationInfo,
    required_object_type = Kratos.Process):

    python_file_name = Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(class_name)
    full_module_name = f"{module_name}.{python_file_name}"

    module = import_module(full_module_name)
    retrieved_object = getattr(module, class_name)(model, parameters, optimization_info)

    # check retrieved object is of the required type
    if not isinstance(retrieved_object, required_object_type):
        raise RuntimeError(f"The retrieved object is of type \"{retrieved_object.__class__.__name__}\" which is not derived from the \"{required_object_type.__name__}\".")
    else:
        return retrieved_object