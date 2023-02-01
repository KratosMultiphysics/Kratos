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

def CallOnAll(list_of_objects: 'list[any]', method: any, *args, **kwargs):
    for obj in list_of_objects:
        getattr(obj, method.__name__)(*args, **kwargs)

def Factory(
    module_name: str,
    class_name: str,
    model: Kratos.Model,
    parameters: Kratos.Parameters,
    optimization_info: OptimizationInfo,
    required_object_type = object):

    python_file_name = ''.join(['_' + c.lower() if c.isupper() else c for c in class_name]).lstrip('_')
    full_module_name = f"{module_name}.{python_file_name}"

    module = import_module(full_module_name)
    retrieved_object = getattr(module, class_name)(model, parameters, optimization_info)

    # check retrieved object is of the required type
    if not isinstance(retrieved_object, required_object_type):
        raise RuntimeError(f"The retrieved object is of type \"{retrieved_object.__class__.__name__}\" which is not of the required type \"{required_object_type.__name__}\".")
    else:
        return retrieved_object