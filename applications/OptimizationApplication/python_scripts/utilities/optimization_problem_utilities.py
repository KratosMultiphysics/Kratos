from importlib import import_module
from typing import Any, Union

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy

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
    for component_type in optimization_problem.GetComponentContainer().keys():
        snake_case_name = Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(component_type.__name__)
        if component_full_name.startswith(snake_case_name):
            return optimization_problem.GetComponent(component_full_name[len(snake_case_name) + 1:], component_type)

    data_container = optimization_problem.GetProblemDataContainer()
    if data_container.HasValue("object") and data_container["object"].HasValue(component_full_name):
        return component_full_name

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

class TensorAdaptorData:
    def __init__(self, tensor_adaptor_path: str, tensor_adaptor: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> None:
        self.tensor_adaptor_path = tensor_adaptor_path
        self.container = tensor_adaptor.GetContainer()

    def GetTensorAdaptorName(self) -> str:
        return self.tensor_adaptor_path[self.tensor_adaptor_path.rfind("/")+1:]

    def GetTensorAdaptorPath(self) -> str:
        return self.tensor_adaptor_path

    def GetContainer(self):
        return self.container

    def GetTensorAdaptor(self, optimization_problem: OptimizationProblem) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        data = optimization_problem.GetProblemDataContainer()[self.tensor_adaptor_path]

        if not isinstance(data, Kratos.TensorAdaptors.DoubleTensorAdaptor):
            raise RuntimeError(f"The data type at \"{self.tensor_adaptor_path}\" changed from {Kratos.TensorAdaptors.DoubleTensorAdaptor.__name__} to {type(data).__class__.__name__}. Found data = {data}")
        if not data.HasContainer():
            raise RuntimeError(f"The data at \"{self.tensor_adaptor_path}\" does not represent a {Kratos.TensorAdaptors.DoubleTensorAdaptor.__name__} with a container. Found data = {data}")
        if data.GetContainer() != self.container:
            raise RuntimeError(f"The container at \"{self.tensor_adaptor_path}\" mismatch with the original container. Found data = {data}")

        return data

class CombinedTensorAdaptorData(TensorAdaptorData):
    def __init__(self, combined_tensor_adaptor_path: str, combined_tensor_adaptor: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor, tensor_adaptor_index: int) -> None:
        super().__init__(combined_tensor_adaptor_path, combined_tensor_adaptor.GetTensorAdaptors()[tensor_adaptor_index])
        self.tensor_adaptor_index = tensor_adaptor_index

    def GetTensorAdaptor(self, optimization_problem: OptimizationProblem) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        data = optimization_problem.GetProblemDataContainer()[self.tensor_adaptor_path]

        if not isinstance(data, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
            raise RuntimeError(f"The data type at \"{self.tensor_adaptor_path}\" changed from {Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor.__name__} to {type(data).__class__.__name__}. Found data = {data}")
        if not data.GetTensorAdaptors()[self.tensor_adaptor_index].HasContainer():
            raise RuntimeError(f"The tensor adaptor at \"{self.tensor_adaptor_path}\" with index = {self.tensor_adaptor_index} does not represent a {Kratos.TensorAdaptors.DoubleTensorAdaptor.__name__} with a container. Found data = {data}")
        if data.GetTensorAdaptors()[self.tensor_adaptor_index].GetContainer() != self.container:
            raise RuntimeError(f"The container from tensor adaptor \"{self.tensor_adaptor_path}\" with index = {self.tensor_adaptor_index} mismatch with the original container. Found data = {data}")

        return data.GetTensorAdaptors()[self.tensor_adaptor_index]

class TensorAdaptorDataManager:
    def __init__(self, optimization_problem: OptimizationProblem):
        self.optimization_problem = optimization_problem
        self.list_of_component_names = ["all"]
        self.list_of_tensor_adaptor_hdf5_outputs = []
        self.create_method = None

    def InitializeHDF5Output(self) -> None:
        # get all the component names at the first writing point
        if len(self.list_of_component_names) == 1 and self.list_of_component_names[0] == "all":
            self.list_of_component_names = GetAllComponentFullNamesWithData(self.optimization_problem)

        list_of_components: 'list[Union[str, ResponseFunction, Control, ExecutionPolicy]]' = []
        for component_name in self.list_of_component_names:
            list_of_components.append(GetComponentHavingDataByFullName(component_name, self.optimization_problem))

        global_values_map = self.optimization_problem.GetProblemDataContainer().GetMap()
        for global_k, global_v in global_values_map.items():
             # first check whether this is part of requested list of components
            found_valid_component = False
            for component in list_of_components:
                component_data = ComponentDataView(component, self.optimization_problem)
                if global_k.startswith(component_data.GetDataPath()):
                     found_valid_component = True
                     break

            # if a valid component is found, add the tensor adaptor
            if found_valid_component:
                if isinstance(global_v, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
                    for i, _ in enumerate(global_v.GetTensorAdaptors()):
                        self._AddContainerTensorAdaptor(CombinedTensorAdaptorData(global_k, global_v, i))
                elif isinstance(global_v, Kratos.TensorAdaptors.DoubleTensorAdaptor) or  \
                     isinstance(global_v, Kratos.TensorAdaptors.IntTensorAdaptor) or \
                     isinstance(global_v, Kratos.TensorAdaptors.BoolTensorAdaptor):
                    self._AddContainerTensorAdaptor(TensorAdaptorData(global_k, global_v))

    def _AddContainerTensorAdaptor(self, tensor_adaptor_data: TensorAdaptorData):
        found_vtu_output = False
        for tensor_adaptor_hdf5_output in self.list_of_tensor_adaptor_hdf5_outputs:
            if tensor_adaptor_hdf5_output.AddTensorAdaptorData(tensor_adaptor_data):
                found_vtu_output = True
                break

        if not found_vtu_output:
            tensor_adaptor_hdf5_output = self.create_method(self.__GetRootModelPart(tensor_adaptor_data.GetContainer()), self.optimization_problem)
            tensor_adaptor_hdf5_output.AddTensorAdaptorData(tensor_adaptor_data)
            self.list_of_tensor_adaptor_hdf5_outputs.append(tensor_adaptor_hdf5_output)
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created tensor adaptor hdf5 output for {tensor_adaptor_hdf5_output.model_part.FullName()}.")

    def __GetRootModelPart(self, container) -> Kratos.ModelPart:
        def get_model_part(container, model_part: Kratos.ModelPart):
            if container in [model_part.Nodes, model_part.Conditions, model_part.Elements]:
                return model_part.GetRootModelPart()

            for sub_model_part_name in model_part.GetSubModelPartNames():
                root_model_part = get_model_part(container, model_part.GetSubModelPart(sub_model_part_name))
                if root_model_part is not None:
                    return root_model_part

            return None

        for model_part_name in self.model.GetModelPartNames():
            root_model_part = get_model_part(container, self.model[model_part_name])
            if root_model_part is not None:
                return root_model_part

        raise RuntimeError(f"No model part contains the provided container.")