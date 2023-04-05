from pathlib import Path
from importlib import import_module
from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
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

ContainerVariableDataHolderUnion = Union[
        KratosOA.HistoricalContainerVariableDataHolder,
        KratosOA.NodalContainerVariableDataHolder,
        KratosOA.ConditionContainerVariableDataHolder,
        KratosOA.ElementContainerVariableDataHolder,
        KratosOA.ConditionPropertiesContainerVariableDataHolder,
        KratosOA.ElementPropertiesContainerVariableDataHolder]

ContainerTypes = Union[
    Kratos.NodesArray,
    Kratos.ConditionsArray,
    Kratos.ElementsArray
]

def ConvertContainers(input_container: ContainerVariableDataHolderUnion, output_container: ContainerVariableDataHolderUnion, number_of_neighbour_entities_for_nodes: ContainerVariableDataHolderUnion):
    if isinstance(input_container, KratosOA.NodalContainerVariableDataHolderBase):
        if isinstance(output_container, KratosOA.NodalContainerVariableDataHolderBase):
            output_container.CopyDataFrom(input_container)
        else:
            KratosOA.ContainerVariableDataHolderUtils.MapNodalVariableDataHolderToContainerVariableDataHolder(output_container, input_container)
    else:
        if isinstance(output_container, KratosOA.NodalContainerVariableDataHolderBase):
            KratosOA.ContainerVariableDataHolderUtils.MapContainerVariableDataHolderToNodalVariableDataHolder(output_container, input_container, number_of_neighbour_entities_for_nodes)
        else:
            raise RuntimeError(f"Unsupported container conversion requested for {input_container} -> {output_container}. Followings are supported: \n\tnodal -> condition\n\tnodal -> element\n\tcondition -> nodal\n\telement -> nodal" )

def WriteCollectiveVariableDataHolderToOptmizationInfo(optimization_info: OptimizationInfo, container: KratosOA.CollectiveVariableDataHolder, pattern: str):
    for sub_container in container.GetVariableDataHolders():
        current_pattern = str(pattern).replace("<model_part_name>", sub_container.GetModelPart().FullName())
        optimization_info.SetValue(current_pattern, sub_container)

def CallOnAll(list_of_objects: 'list[any]', method: any, *args, **kwargs):
    for obj in list_of_objects:
        getattr(obj, method.__name__)(*args, **kwargs)

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
