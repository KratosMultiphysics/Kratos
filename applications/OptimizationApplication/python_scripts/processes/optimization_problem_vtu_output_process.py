from typing import Any, Union
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Any:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemVtuOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemVtuOutputProcess(model, parameters["settings"], optimization_problem)

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

class TensorAdaptorVtuOutput:
    def __init__(self, parameters: 'dict[str, Any]', model_part: Kratos.ModelPart, optimization_problem: OptimizationProblem):
        self.model_part = model_part
        self.optimization_problem = optimization_problem
        self.output_file_name_prefix: str = parameters["output_file_name_prefix"]

        if parameters["save_output_files_in_folder"]:
            self.output_path = Path(parameters["output_path"])
            if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
            # now create the output path
            Kratos.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path))
        else:
            self.output_path = Path(".")

        self.vtu_output: Kratos.Future.VtuOutput = Kratos.Future.VtuOutput(model_part, parameters["configuration"], parameters["writer_format"], parameters["precision"], output_sub_model_parts=True, echo_level=parameters["echo_level"], write_ids=parameters["write_ids"])
        self.list_of_tensor_adaptor_data: 'list[TensorAdaptorData]' = []

    def AddTensorAdaptorData(self, tensor_adaptor_data: TensorAdaptorData) -> bool:
        if tensor_adaptor_data.GetContainer() in self.vtu_output.GetOutputContainerList():
            self.list_of_tensor_adaptor_data.append(tensor_adaptor_data)
            return True
        return False

    def WriteOutput(self):
        for tensor_adaptor_data in self.list_of_tensor_adaptor_data:
            self.vtu_output.EmplaceTensorAdaptor(tensor_adaptor_data.GetTensorAdaptorName(), tensor_adaptor_data.GetTensorAdaptor(self.optimization_problem))

        output_file_name = self.output_file_name_prefix
        output_file_name = output_file_name.replace("<model_part_full_name>", self.model_part.FullName())
        output_file_name = output_file_name.replace("<model_part_name>", self.model_part.Name)
        self.vtu_output.PrintOutput(str(self.output_path / output_file_name))

class OptimizationProblemVtuOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "file_name"                   : "<model_part_full_name>",
                "file_format"                 : "binary",
                "output_path"                 : "Optimization_Results",
                "save_output_files_in_folder" : true,
                "write_ids"                   : false,
                "write_deformed_configuration": false,
                "list_of_output_components"   : ["all"],
                "output_precision"            : 7,
                "output_interval"             : 1,
                "echo_level"                  : 0
            }
            """
        )

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__()

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.optimization_problem = optimization_problem

        self.echo_level = parameters["echo_level"].GetInt()
        self.file_name = parameters["file_name"].GetString()
        self.output_path = parameters["output_path"].GetString()
        self.save_output_files_in_folder = parameters["save_output_files_in_folder"].GetBool()
        self.output_precision = parameters["output_precision"].GetInt()
        self.write_ids = parameters["write_ids"].GetBool()
        file_format = parameters["file_format"].GetString()
        if file_format == "ascii":
            self.writer_format = Kratos.Future.VtuOutput.ASCII
        elif file_format == "binary":
            self.writer_format = Kratos.Future.VtuOutput.BINARY
        elif file_format == "raw":
            self.writer_format = Kratos.Future.VtuOutput.RAW
        elif file_format == "compressed_raw":
            self.writer_format = Kratos.Future.VtuOutput.COMPRESSED_RAW
        else:
            raise RuntimeError(f"Only supports \"ascii\", \"binary\", \"raw\", and \"compressed_raw\" file_format. [ provided file_format = \"{file_format}\" ].")

        if parameters["write_deformed_configuration"].GetBool():
            self.configuration = Kratos.Configuration.Current
        else:
            self.configuration = Kratos.Configuration.Initial

        self.list_of_component_names = parameters["list_of_output_components"].GetStringArray()
        self.list_of_tensor_adaptor_vtu_outputs: 'list[TensorAdaptorVtuOutput]' = []
        self.initialized_vtu_outputs = False

    def PrintOutput(self) -> None:
        if not self.initialized_vtu_outputs:
            self.InitializeVtuOutputIO()
            self.initialized_vtu_outputs = True

        for tensor_adaptor_vtu_output in self.list_of_tensor_adaptor_vtu_outputs:
            tensor_adaptor_vtu_output.WriteOutput()

    def InitializeVtuOutputIO(self) -> None:
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
        for tensor_adaptor_vtu_output in self.list_of_tensor_adaptor_vtu_outputs:
            if tensor_adaptor_vtu_output.AddTensorAdaptorData(tensor_adaptor_data):
                found_vtu_output = True
                break

        if not found_vtu_output:
            vtu_parameters = {
                "output_file_name_prefix": self.file_name,
                "configuration": self.configuration,
                "writer_format": self.writer_format,
                "precision": self.output_precision,
                "save_output_files_in_folder": self.save_output_files_in_folder,
                "output_path": self.output_path,
                "echo_level": self.echo_level,
                "write_ids": self.write_ids
            }

            tensor_adaptor_vtu_output = TensorAdaptorVtuOutput(vtu_parameters, self.__GetRootModelPart(tensor_adaptor_data.GetContainer()), self.optimization_problem)
            tensor_adaptor_vtu_output.AddTensorAdaptorData(tensor_adaptor_data)
            self.list_of_tensor_adaptor_vtu_outputs.append(tensor_adaptor_vtu_output)
            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created tensor adaptor vtu output for {tensor_adaptor_vtu_output.vtu_output.GetModelPart().FullName()}.")

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