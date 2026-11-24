from pathlib import Path
from typing import Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.kratos_utilities import GetListOfAvailableApplications
from KratosMultiphysics.kratos_utilities import GetKratosMultiphysicsPath

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
