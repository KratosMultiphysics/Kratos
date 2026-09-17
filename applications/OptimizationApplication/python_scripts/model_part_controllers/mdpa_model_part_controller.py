import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"MdpaModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return MdpaModelPartController(model, parameters["settings"])

class MdpaModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "model_part_name": "",
            "input_filename" : "",
            "domain_size"    : -1,
            "read_data"      : false
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
            raise RuntimeError("\"domain_size\"  should be either 1, 2 or 3." + str(parameters))

        self.model_part = model.CreateModelPart(model_part_name)
        self.read_data = parameters["read_data"].GetBool()

    def SupportsRestart(self) -> bool:
        return True

    def ImportModelPart(self) -> None:
        if self.model_part.NumberOfNodes() > 0:
            # OptimizationAnalysis._CreateRestart() already loaded this model part directly from
            # the restart checkpoint, synchronously, before Initialize() ever calls this method --
            # skip the normal mdpa read.
            # Mirrors Kratos.RestartUtility.LoadRestart(): downstream solvers (e.g.
            # MechanicalSolver.PrepareModelPart()) gate re-reading materials/constitutive laws on
            # IS_RESTARTED -- skipping it is required, not just an optimization, since the
            # restored model part's Properties/constitutive laws would otherwise be clobbered by
            # a second, independent materials-import on top of the deserialized ones.
            self.model_part.ProcessInfo[Kratos.IS_RESTARTED] = True
            self.model_part.ProcessInfo[Kratos.LOAD_RESTART] = self.model_part.ProcessInfo[Kratos.STEP] + 1
            return

        if self.read_data:
            Kratos.ModelPartIO(self.input_filename, Kratos.ModelPartIO.READ).ReadModelPart(self.model_part)
        else:
            Kratos.ModelPartIO(self.input_filename, Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(self.model_part)

        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = self.domain_size

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part
