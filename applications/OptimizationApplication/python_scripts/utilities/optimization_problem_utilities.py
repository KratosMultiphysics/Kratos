from pathlib import Path
from importlib import import_module
from typing import Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def OptimizationComponentFactory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    if not parameters.Has("type"):
        raise RuntimeError(f"Components created from OptimizationComponentFactory require the \"type\" [ provided parameters = {parameters}].")

    python_type = parameters["type"].GetString()

    if not parameters.Has("module") or parameters["module"].GetString() == "":
        # in the case python type comes without a module
        # as in the case python_type is in the sys path or the current working directory.
        full_module_name = python_type
    else:
        # in the case python type comes with a module.
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

def GetComponentValueByFullName(component_data_view: ComponentDataView, value_name: str, step_index = 0) -> Any:
    value_info = value_name.split(":")
    if len(value_info) != 2 or value_info[1] not in ["historical", "non_historical"]:
        raise RuntimeError(f"The value_name should contain the name of the value and either it is located at buffered data or un buffered data. Ex: \"objective_Value:historical\" for buffered container or \"reference_value:non_historical\" for unbuffered container [ value_name = \"{value_name}\" ].")

    if value_info[1] == "historical":
        data = component_data_view.GetBufferedData()
    elif value_info[1] == "non_historical":
        data = component_data_view.GetUnBufferedData()

    if not data.HasValue(value_info[0], step_index):
        raise RuntimeError(f"The value = \"{value_info[0]}\" is not found in the component = \"{component_data_view.GetComponentName()}\" [ value_name = \"{value_name}\", step_index = {step_index}, data = {data}, component_data_view = {component_data_view}")

    return data.GetValue(value_info[0], step_index)

def SetComponentValueByFullName(component_data_view: ComponentDataView, value_name: str, value: Any, step_index = 0, overwrite = False) -> None:
    value_info = value_name.split(":")
    if len(value_info) != 2 or value_info[1] not in ["historical", "non_historical"]:
        raise RuntimeError(f"The value_name should contain the name of the value and either it is located at buffered data or un buffered data. Ex: \"objective_Value:historical\" for buffered container or \"reference_value:non_historical\" for unbuffered container [ value_name = \"{value_name}\" ].")

    if value_info[1] == "historical":
        data = component_data_view.GetBufferedData()
    elif value_info[1] == "non_historical":
        data = component_data_view.GetUnBufferedData()

    data.SetValue(value_info[0], value, step_index, overwrite)

def OutputGradientFields(response: ResponseRoutine, optimization_problem: OptimizationProblem, is_gradients_computed: bool) -> None:
    unbuffered_data = ComponentDataView(response.GetReponse(), optimization_problem).GetUnBufferedData()
    if is_gradients_computed:
        # save the physical gradients for post processing in unbuffered data container.
        for physical_var, physical_gradient in response.GetRequiredPhysicalGradients().items():
            variable_name = f"d{response.GetResponseName()}_d{physical_var.Name()}"
            for physical_gradient_ta in physical_gradient.GetTensorAdaptors():
                unbuffered_data.SetValue(variable_name, Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_gradient_ta), overwrite=True)

        # save the filtered gradients for post processing in unbuffered data container.
        for gradient_container_ta, control in zip(response.GetMappedGradients().GetTensorAdaptors(), response.GetMasterControl().GetListOfControls()):
            variable_name = f"d{response.GetResponseName()}_d{control.GetName()}"
            unbuffered_data.SetValue(variable_name, Kratos.TensorAdaptors.DoubleTensorAdaptor(gradient_container_ta), overwrite=True)
    else:
        # save the physical gradients for post processing in unbuffered data container.
        for physical_var, physical_gradient in response.GetRequiredPhysicalGradients().items():
            variable_name = f"d{response.GetResponseName()}_d{physical_var.Name()}"
            for physical_gradient_ta in physical_gradient.GetTensorAdaptors():
                temp_ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_gradient_ta)
                temp_ta.data[:] = 0.0
                unbuffered_data.SetValue(variable_name, temp_ta, overwrite=True)

        # save the filtered gradients for post processing in unbuffered data container.
        for control in response.GetMasterControl().GetListOfControls():
            variable_name = f"d{response.GetResponseName()}_d{control.GetName()}"
            unbuffered_data.SetValue(variable_name, control.GetEmptyField(), overwrite=True)
