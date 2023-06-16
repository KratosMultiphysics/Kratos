from typing import Any
from enum import Enum
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class ModelPartUtilities:
    class OperationType(Enum):
        UNION = 1,
        INTERSECT = 2,
        SUBSTRACT = 3

    @staticmethod
    def __GenerateUniqueModelPartName(prefix: str, list_of_names: 'list[str]', add_neighbours: bool) -> str:
        if not all([name.find("#") == -1 for name in list_of_names]):
            raise RuntimeError(f"The provided model part names has \"#\" which is invalid. Please remove all occurances of \"#\" character. Provided model part names:\n\t" + "\n\t".join(list_of_names))

        sorted_names = sorted(list_of_names)
        post_fix = "IN" if add_neighbours else "EN"
        return (f"{prefix}_{post_fix}_" + "_".join(sorted_names)).replace(".", "-")

    @staticmethod
    def GetOperatingModelPart(operation_type: OperationType, main_model_part: Kratos.ModelPart, operation_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> Kratos.ModelPart:
        model_part_names = [model_part.FullName() for model_part in operation_model_parts]
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName(operation_type.name, model_part_names, add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return main_model_part.GetSubModelPart(unique_name)
        else:
            unique_model_part = main_model_part.CreateSubModelPart(unique_name)
            KratosOA.ModelPartUtils.LogModelPartStatus(unique_model_part, f"ModelPartUtilities_created#{operation_type.name}#{add_neighbours:d}#" + "#".join(model_part_names))
            return unique_model_part

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
                        elif operation_type == "SUBSTRACT":
                            Kratos.ModelPartOperationUtilities.Substract(model_part, main_model_part, operation_model_parts, add_neighbours)
                        else:
                            raise RuntimeError(f"Unsupported operation type found [ operation_type = {operation_type}, model_part = {model_part.FullName()}]. Followings are list of statuses:\n\t" + "\n\t".join(statuses))

                        KratosOA.ModelPartUtils.LogModelPartStatus(model_part, filled_status)
                        return True
                    else:
                        raise RuntimeError(f"The model part {model_part.FullName()} cannot be properly constructed. Followings are list of statuses:\n\t" + "\n\t".join(statuses))

        return False

    @staticmethod
    def UnionModelParts(main_model_part: Kratos.ModelPart, union_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> 'tuple[bool, Kratos.ModelPart]':
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName("UNION", [model_part.FullName() for model_part in union_model_parts], add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return False, main_model_part.GetSubModelPart(unique_name)
        else:
            unique_model_part = main_model_part.CreateSubModelPart(unique_name)
            Kratos.ModelPartOperationUtilities.Union(unique_model_part, main_model_part, union_model_parts, add_neighbours)
            return True, unique_model_part

    @staticmethod
    def IntersectModelParts(main_model_part: Kratos.ModelPart, intersecting_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> 'tuple[bool, Kratos.ModelPart]':
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName("INTERSECT", [model_part.FullName() for model_part in intersecting_model_parts], add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return False, main_model_part.GetSubModelPart(unique_name)
        else:
            unique_model_part = main_model_part.CreateSubModelPart(unique_name)
            Kratos.ModelPartOperationUtilities.Intersect(unique_model_part, main_model_part, intersecting_model_parts, add_neighbours)
            return True, unique_model_part

    @staticmethod
    def SubstractModelParts(main_model_part: Kratos.ModelPart, substracting_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> 'tuple[bool, Kratos.ModelPart]':
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName("SUBSTRACT", [model_part.FullName() for model_part in substracting_model_parts], add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return False, main_model_part.GetSubModelPart(unique_name)
        else:
            unique_model_part = main_model_part.CreateSubModelPart(unique_name)
            Kratos.ModelPartOperationUtilities.Substract(unique_model_part, main_model_part, substracting_model_parts, add_neighbours)
            return True, unique_model_part

    @staticmethod
    def GetMergedMap(main_model_part: Kratos.ModelPart, input_dict: 'dict[Any, KratosOA.CollectiveExpression]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            merging_model_parts = [container_expression.GetModelPart() for container_expression in v.GetContainerExpressions()]
            merged_model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, main_model_part, merging_model_parts, add_neghbours)
            ModelPartUtilities.ExecuteOperationOnModelPart(merged_model_part)
            result[k] = merged_model_part

        return result

    @staticmethod
    def GetIntersectedMap(main_model_part: Kratos.ModelPart, input_dict: 'dict[Any, Kratos.ModelPart]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            intersected_model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.INTERSECT, main_model_part, [v], add_neghbours)
            ModelPartUtilities.ExecuteOperationOnModelPart(intersected_model_part)
            result[k] = intersected_model_part

        return result