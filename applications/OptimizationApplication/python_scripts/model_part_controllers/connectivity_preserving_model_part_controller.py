import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"ConnectivityPreservingModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return ConnectivityPreservingModelPartController(model, parameters["settings"], optimization_problem)

class ConnectivityPreservingModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "transformation_settings": [
                {
                    "source_model_part_name"     : "",
                    "destination_model_part_name": "",
                    "destination_element_name"   : "",
                    "destination_condition_name" : ""
                }
            ]
        }""")

        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.optimization_problem = optimization_problem

        for sub_model_part_settings in self.parameters["transformation_settings"].values():
            sub_model_part_settings.ValidateAndAssignDefaults(default_settings["transformation_settings"].values()[0])

    def ImportModelPart(self) -> None:
        if self.optimization_problem.IsRestarted():
            # the destination model parts come back populated from the whole-Model restore;
            # re-running ConnectivityPreserveModeler on them would duplicate entities.
            return

        connectivity_preserve_modeller = Kratos.ConnectivityPreserveModeler()
        for transformation_settings in self.parameters["transformation_settings"].values():
            # get source and destination model parts
            source_model_part = self.model[transformation_settings["source_model_part_name"].GetString()]
            destination_model_part_name: str = transformation_settings["destination_model_part_name"].GetString()
            if not self.model.HasModelPart(destination_model_part_name):
                # if the destination model part is not there, then check whether the destination root model part is existing
                destination_root_model_part_name = destination_model_part_name.split(".")[0]
                if not self.model.HasModelPart(destination_root_model_part_name):
                    # destination root model part is also not there. so we can safely assign the source model part solution variables list
                    destination_root_model_part = self.model.CreateModelPart(destination_root_model_part_name)
                else:
                    # if it is existing, then check whether it has a subset of the variables list.
                    destination_root_model_part = self.model[destination_root_model_part_name]
                    if not KratosOA.OptimizationUtils.IsSolutionStepVariablesListASubSet(source_model_part, destination_root_model_part):
                        msg = f"The destination root model part solution step variables at {destination_root_model_part.FullName()} is not a sub set of {source_model_part.FullName()}"
                        msg += f"\n\t{destination_root_model_part.FullName()} model part variables list:\n\t\t" + "\n\t\t".join(KratosOA.OptimizationUtils.GetSolutionStepVariableNamesList(destination_root_model_part))
                        msg += f"\n\t{source_model_part.FullName()} model part variables list:\n\t\t" + "\n\t\t".join(KratosOA.OptimizationUtils.GetSolutionStepVariableNamesList(source_model_part))
                        raise RuntimeError(msg)
                KratosOA.OptimizationUtils.SetSolutionStepVariablesList(destination_root_model_part, source_model_part)
                destination_model_part = self.model.CreateModelPart(destination_model_part_name)
            else:
                destination_model_part = self.model[destination_model_part_name]

            element_name = transformation_settings["destination_element_name"].GetString()
            condition_name = transformation_settings["destination_condition_name"].GetString()

            if condition_name != "" and element_name != "":
                connectivity_preserve_modeller.GenerateModelPart(source_model_part, destination_model_part, element_name, condition_name)
            elif element_name != "":
                connectivity_preserve_modeller.GenerateModelPart(source_model_part, destination_model_part, element_name)
            elif condition_name != "":
                connectivity_preserve_modeller.GenerateModelPart(source_model_part, destination_model_part, condition_name)
            else:
                raise RuntimeError("The element_name and the condition name cannot be empty.")

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Duplicated {source_model_part.FullName()} in to {destination_model_part.FullName()}.")

    def GetModelPart(self) -> Kratos.ModelPart:
        return None