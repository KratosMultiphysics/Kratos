import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"ExternalResponseFunctionModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return ExternalResponseFunctionModelPartController(model, parameters["settings"])

class ExternalResponseFunctionModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "model_part_name"        : "",
            "number_of_design_points": 0
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        model_part_name = parameters["model_part_name"].GetString()
        if model_part_name == "":
            raise RuntimeError("Empty \"model_part_name\" is not allowed which is given with following parameters:\n" + str(parameters))

        self.model_part = model.CreateModelPart(model_part_name)
        self.number_of_nodes = parameters["number_of_design_points"].GetInt()

    def ImportModelPart(self) -> None:
        for i in range(self.number_of_nodes):
            self.model_part.CreateNewNode(i + 1, i + 1, 0, 0)

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part
