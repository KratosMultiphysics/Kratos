from importlib import import_module
from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

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
    optimization_info: OptimizationInfo = None,
    required_object_type = Kratos.Process):

    python_file_name = ''.join(['_' + c.lower() if c.isupper() else c for c in class_name]).lstrip('_')
    full_module_name = f"{module_name}.{python_file_name}"

    module = import_module(full_module_name)
    if optimization_info is not None:
        retrieved_object = getattr(module, class_name)(model, parameters, optimization_info)
    else:
        retrieved_object = getattr(module, class_name)(model, parameters)

    # check retrieved object is of the required type
    if not isinstance(retrieved_object, required_object_type):
        raise RuntimeError(f"The retrieved object is of type \"{retrieved_object.__class__.__name__}\" which is not derrived from the \"{required_object_type.__name__}\".")
    else:
        return retrieved_object
