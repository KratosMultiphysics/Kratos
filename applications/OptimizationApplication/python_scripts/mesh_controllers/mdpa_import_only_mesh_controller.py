import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.mesh_controllers.mesh_controller import MeshController
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class MdpaImportOnlyMeshController(MeshController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        default_settings = Kratos.Parameters("""{
            "model_part_name": "",
            "input_filename" : "",
            "domain_size"    : -1
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        model_part_name = parameters["model_part_name"].GetString()
        if model_part_name == "":
            raise RuntimeError("Empty \"model_part_name\" is not allowed which is given with following parameters:\n" + str(parameters))

        self.input_filename = parameters["input_filename"].GetString()
        if self.input_filename == "":
            raise RuntimeError("Empty \"input_filename\" is not allowed which is given with following parameters:\n" + str(parameters))

        self.domain_size = parameters["domain_size"].GetInt()
        if self.domain_size <= 0 or self.domain_size > 3:
            raise RuntimeError("\"domain_size\"  should be in [1, 3] range." + str(parameters))

        self.model_part = self.model.CreateModelPart(model_part_name)

    def ImportModelPart(self):
        Kratos.ModelPartIO(self.input_filename, Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(self.model_part)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = self.domain_size