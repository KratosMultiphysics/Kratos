import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import TensorAdaptorData
from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_field_output import TensorAdaptorOutput, OptimizationProblemFieldOutput

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemHDF5OutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemHDF5OutputProcess(model, parameters["settings"], optimization_problem)

class TensorAdaptorHDF5Output(TensorAdaptorOutput):
    def __init__(self, model_part: Kratos.ModelPart, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):

        default_parameters = Kratos.Parameters("""{
            "file_name"       : "PLEASE_SPECIFY_HDF5_FILENAME",
            "file_access_mode": "exclusive",
            "echo_level"      : 0
        }""")

        self.model_part = model_part
        self.optimization_problem = optimization_problem
        self.output_parameters = parameters["output_file_settings"]
        self.output_path_prefix = parameters["hdf5_output_path_prefix"].GetString()
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

class OptimizationProblemHDF5OutputProcess(OptimizationProblemFieldOutput):
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

    def _CreateTensorAdaptorOutput(self, tensor_adaptor_data: TensorAdaptorData) -> TensorAdaptorOutput:
        return TensorAdaptorHDF5Output(self._GetModelPart(tensor_adaptor_data.GetContainer()), self.parameters.Clone(), self.optimization_problem)
