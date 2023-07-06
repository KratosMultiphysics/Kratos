from typing import Any
from enum import Enum
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class ModelPartUtilities:
    class OperationType(Enum):
        UNION = 1,
        INTERSECT = 2

    @staticmethod
    def __GenerateUniqueIdentifier(prefix: str, list_of_names: 'list[str]', add_neighbours: bool) -> str:
        if not all([name.find("#") == -1 for name in list_of_names]):
            raise RuntimeError(f"The provided model part names has \"#\" which is invalid. Please remove all occurances of \"#\" character. Provided model part names:\n\t" + "\n\t".join(list_of_names))

        sorted_names = sorted(list_of_names)
        post_fix = "IN" if add_neighbours else "EN"
        return (f"{prefix}_{post_fix}_" + "_".join(sorted_names)).replace(".", "-")

    @staticmethod
    def GetOperatingModelPart(operation_type: OperationType, suggested_model_part_name: str, operation_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> Kratos.ModelPart:
        if len(operation_model_parts) == 0:
            raise RuntimeError("No operating model parts are provided.")

        # check whether is there a need to perform binary operations
        if operation_type == ModelPartUtilities.OperationType.UNION or operation_type == ModelPartUtilities.OperationType.INTERSECT:
            if len(set(operation_model_parts)) == 1:
                # all the operation model parts are the same, hence no need of doing an operation.
                return operation_model_parts[0]

        model_part_names = [model_part.FullName() for model_part in operation_model_parts]
        model_part_names = sorted(model_part_names)
        root_model_part: Kratos.ModelPart = operation_model_parts[0].GetModel()[model_part_names[0]].GetRootModelPart()
        operation_identifier = f"ModelPartUtilities_created#{operation_type.name}#{add_neighbours:d}#" + "#".join(model_part_names)

        for sub_model_part in root_model_part.SubModelParts:
            if KratosOA.ModelPartUtils.CheckModelPartStatus(sub_model_part, operation_identifier):
                Kratos.Logger.PrintInfo("ModelPartUtilities", f"Using {sub_model_part.FullName()} with the same operation instead of suggested model part with name = \"{suggested_model_part_name}\".")
                return sub_model_part

        if root_model_part.HasSubModelPart(suggested_model_part_name):
            # this means, it already has a model part with suggeted name, but
            # it does not match the operation identifier. So throw an error
            raise RuntimeError(f"Found an already existing submodel part named \"{suggested_model_part_name}\" in {root_model_part.FullName()} without the required operation identifier = \"{operation_identifier}\".")

        # if the operation identifier is not found with multiple operation_model_partss, then have to create it
        # the sub model part
        sub_model_part = root_model_part.CreateSubModelPart(suggested_model_part_name)
        KratosOA.ModelPartUtils.LogModelPartStatus(sub_model_part, operation_identifier)
        Kratos.Logger.PrintInfo("ModelPartUtilities", f"Created sub model part \"{sub_model_part.FullName()}\".")
        return sub_model_part

    @staticmethod
    def ExecuteOperationOnModelPart(model_part: Kratos.ModelPart) -> bool:
        statuses = KratosOA.ModelPartUtils.GetModelPartStatusLog(model_part)
        for status in statuses:
            status_data = status.split("#")
            if len(status_data) > 0 and status_data[0].startswith("ModelPartUtilities_created"):
                # found a model part which was created using this utility. check whether it is
                # filled with necessary entitites.

                filled_status = "ModelPartUtilities_filled#" + "#".join(status_data[1:])
                if not KratosOA.ModelPartUtils.CheckModelPartStatus(model_part, filled_status):
                    # the model part is not filled with the corresponding entities
                    # hence, fill it now
                    if len(status_data) > 3:
                        model = model_part.GetModel()
                        operation_type = status_data[1]
                        add_neighbours = bool(status_data[2])
                        operation_model_parts = [model[name] for name in status_data[3:]]
                        main_model_part = model_part.GetParentModelPart()

                        if operation_type == "UNION":
                            Kratos.ModelPartOperationUtilities.Union(model_part, main_model_part, operation_model_parts, add_neighbours)
                        elif operation_type == "INTERSECT":
                            Kratos.ModelPartOperationUtilities.Intersect(model_part, main_model_part, operation_model_parts, add_neighbours)
                        else:
                            raise RuntimeError(f"Unsupported operation type found [ operation_type = {operation_type}, model_part = {model_part.FullName()}]. Followings are list of statuses:\n\t" + "\n\t".join(statuses))

                        KratosOA.ModelPartUtils.LogModelPartStatus(model_part, filled_status)
                        return True
                    else:
                        raise RuntimeError(f"The model part {model_part.FullName()} cannot be properly constructed. Followings are list of statuses:\n\t" + "\n\t".join(statuses))

        return False

    @staticmethod
    def GetMergedMap(input_dict: 'dict[Any, KratosOA.CollectiveExpression]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            merging_model_parts = [container_expression.GetModelPart() for container_expression in v.GetContainerExpressions()]
            uniqe_identifier_name = ModelPartUtilities.__GenerateUniqueIdentifier("UNION", [model_part.FullName() for model_part in merging_model_parts], add_neghbours)
            merged_model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, uniqe_identifier_name, merging_model_parts, add_neghbours)
            ModelPartUtilities.ExecuteOperationOnModelPart(merged_model_part)
            result[k] = merged_model_part

        return result

    @staticmethod
    def GetIntersectedMap(main_model_part: Kratos.ModelPart, input_dict: 'dict[Any, Kratos.ModelPart]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            uniqe_identifier_name = ModelPartUtilities.__GenerateUniqueIdentifier("INTERSECT", [main_model_part.FullName(), v.FullName()], add_neghbours)
            intersected_model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.INTERSECT, uniqe_identifier_name, [main_model_part, v], add_neghbours)
            ModelPartUtilities.ExecuteOperationOnModelPart(intersected_model_part)
            result[k] = intersected_model_part

        return result