from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import TensorAdaptorData
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_field_output_process import TensorAdaptorOutput, OptimizationProblemFieldOutputProcess

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemVtuOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemVtuOutputProcess(model, parameters["settings"], optimization_problem)

class TensorAdaptorVtuOutput(TensorAdaptorOutput):
    def __init__(self,  model_part: Kratos.ModelPart, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model_part = model_part.GetRootModelPart()
        self.optimization_problem = optimization_problem

        if parameters["save_output_files_in_folder"].GetBool():
            self.output_path = Path(parameters["output_path"].GetString())
            if not self.model_part.ProcessInfo[Kratos.IS_RESTARTED]:
                kratos_utils.DeleteDirectoryIfExisting(str(self.output_path))
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()
            # now create the output path
            Kratos.FilesystemExtensions.MPISafeCreateDirectories(str(self.output_path))
        else:
            self.output_path = Path(".")

        file_format = parameters["file_format"].GetString()
        if file_format == "ascii":
            self.writer_format = Kratos.VtuOutput.ASCII
        elif file_format == "binary":
            self.writer_format = Kratos.VtuOutput.BINARY
        elif file_format == "raw":
            self.writer_format = Kratos.VtuOutput.RAW
        elif file_format == "compressed_raw":
            self.writer_format = Kratos.VtuOutput.COMPRESSED_RAW
        else:
            raise RuntimeError(f"Only supports \"ascii\", \"binary\", \"raw\", and \"compressed_raw\" file_format. [ provided file_format = \"{file_format}\" ].")

        self.output_file_name_prefix = parameters["file_name"].GetString()
        self.vtu_output: Kratos.VtuOutput = Kratos.VtuOutput(
                                                model_part,
                                                not parameters["write_deformed_configuration"].GetBool(),
                                                self.writer_format,
                                                parameters["output_precision"].GetInt(),
                                                output_sub_model_parts=True,
                                                echo_level=parameters["echo_level"].GetInt(),
                                                write_ids=parameters["write_ids"].GetBool())

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

class OptimizationProblemVtuOutputProcess(OptimizationProblemFieldOutputProcess):
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

    def _CreateTensorAdaptorOutput(self, tensor_adaptor_data: TensorAdaptorData) -> TensorAdaptorOutput:
        return TensorAdaptorVtuOutput(self._GetModelPart(tensor_adaptor_data.GetContainer()), self.parameters.Clone(), self.optimization_problem)
