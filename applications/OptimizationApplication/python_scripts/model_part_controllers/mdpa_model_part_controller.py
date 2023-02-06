import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class MdpaModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "model_part_name": "",
            "input_filename" : "",
            "domain_size"    : -1,
            "echo_level"     : 0
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        model_part_name = parameters["model_part_name"].GetString()
        if model_part_name == "":
            raise RuntimeError("Empty \"model_part_name\" is not allowed which is given with following parameters:\n" + str(parameters))

        self.input_filename = parameters["input_filename"].GetString()
        if self.input_filename == "":
            raise RuntimeError("Empty \"input_filename\" is not allowed which is given with following parameters:\n" + str(parameters))

        self.domain_size = parameters["domain_size"].GetInt()
        if self.domain_size not in [1, 2, 3]:
            raise RuntimeError("\"domain_size\"  should be in either 1, 2 or 3." + str(parameters))

        self.model_part = model.CreateModelPart(model_part_name)
        self.optimization_info = optimization_info
        self.echo_level = parameters["echo_level"].GetInt()

    def ImportModelPart(self):
        Kratos.ModelPartIO(self.input_filename, Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(self.model_part)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = self.domain_size

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Read model part {self.model_part.FullName()} from {self.input_filename}.")

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def InitializeSolutionStep(self):
        self.model_part.ProcessInfo[Kratos.STEP] = self.optimization_info["step"]
        self.model_part.ProcessInfo[Kratos.TIME] = self.optimization_info["step"]