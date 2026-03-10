from typing import Any, Union
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import TensorAdaptorData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import CombinedTensorAdaptorData

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemHDF5OutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemHDF5OutputProcess(model, parameters["settings"], optimization_problem)

class TensorAdaptorHDF5Output:
    def __init__(self, output_path_prefix: str, parameters: Kratos.Parameters, model_part: Kratos.ModelPart, optimization_problem: OptimizationProblem):
        default_parameters = Kratos.Parameters("""{
            "file_name"       : "PLEASE_SPECIFY_HDF5_FILENAME",
            "file_access_mode": "exclusive",
            "echo_level"      : 0
        }""")

        self.model_part = model_part
        self.optimization_problem = optimization_problem
        self.output_parameters = parameters
        self.output_path_prefix = output_path_prefix
        self.list_of_tensor_adaptor_data: 'list[TensorAdaptorData]' = []

        self.output_parameters.ValidateAndAssignDefaults(default_parameters)
        self.output_parameters["file_name"].SetString(self.output_parameters["file_name"].GetString().replace("<model_part_full_name>", self.model_part.FullName()))
        self.output_parameters["file_name"].SetString(self.output_parameters["file_name"].GetString().replace("<model_part_name>", self.model_part.Name))
        self.output_path_prefix = self.output_path_prefix.replace("<model_part_full_name>", self.model_part.FullName())
        self.output_path_prefix = self.output_path_prefix.replace("<model_part_name>", self.model_part.Name)

    def AddTensorAdaptorData(self, tensor_adaptor_data: TensorAdaptorData) -> bool:
        if tensor_adaptor_data.GetContainer() in [self.model_part.Nodes, self.model_part.Conditions, self.model_part.Elements]:
            self.list_of_tensor_adaptor_data.append(tensor_adaptor_data)
            return True
        return False

    def WriteOutput(self):
        current_output_parameters = self.output_parameters.Clone()
        current_output_parameters["file_name"].SetString(current_output_parameters["file_name"].GetString().replace("<step>", str(self.optimization_problem.GetStep())))

        tensor_io_settings = Kratos.Parameters("""{ "prefix": "" }""")
        tensor_io_settings["prefix"].SetString(self.output_path_prefix.replace("<step>", str(self.optimization_problem.GetStep())))

        with OpenHDF5File(current_output_parameters, self.model_part) as hdf5_file:
            tensor_io = KratosHDF5.TensorAdaptorIO(tensor_io_settings, hdf5_file)
            for tensor_adaptor_data in self.list_of_tensor_adaptor_data:
                tensor_io.Write(tensor_adaptor_data.GetTensorAdaptorName(), tensor_adaptor_data.GetTensorAdaptor(self.optimization_problem))

class OptimizationProblemHDF5OutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "list_of_output_components": ["all"],
                "hdf5_output_path_prefix"  : "Optimization_Results",
                "output_interval"          : 1,
                "echo_level"               : 0,
                "output_file_settings"     : {
                    "file_name"       : "PLEASE_SPECIFY_HDF5_FILENAME",
                    "file_access_mode": "exclusive",
                    "echo_level"      : 0
                }
            }
            """
        )

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__()

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.optimization_problem = optimization_problem
        self.list_of_tensor_adaptor_hdf5_outputs: 'list[TensorAdaptorHDF5Output]' = []
        self.hdf5_output_path_prefix = ""
        self.echo_level = 1
        self.output_file_settings = Kratos.Parameters()

    def PrintOutput(self) -> None:
        if not self.initialized_hdf5_outputs:
            self.InitializeHDF5Output()
            self.initialized_hdf5_outputs = True

        for tensor_adaptor_hdf5_output in self.list_of_tensor_adaptor_hdf5_outputs:
            tensor_adaptor_hdf5_output.WriteOutput()

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
            tensor_adaptor_hdf5_output = TensorAdaptorHDF5Output(self.hdf5_output_path_prefix, self.output_file_settings.Clone(), self.__GetRootModelPart(tensor_adaptor_data.GetContainer()), self.optimization_problem)
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