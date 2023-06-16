from typing import Any
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class ModelPartUtilities:
    @staticmethod
    def __GenerateUniqueModelPartName(prefix: str, list_of_names: 'list[str]', add_neighbours: bool) -> str:
        sorted_names = sorted(list_of_names)
        post_fix = "IN" if add_neighbours else "EN"
        return (f"{prefix}_{post_fix}_" + "_".join(sorted_names)).replace(".", "-")

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
            _, merged_model_part = ModelPartUtilities.UnionModelParts(main_model_part, merging_model_parts, add_neghbours)
            result[k] = merged_model_part

        return result

    @staticmethod
    def GetIntersectedMap(main_model_part: Kratos.ModelPart, input_dict: 'dict[Any, Kratos.ModelPart]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
        result: 'dict[Any, Kratos.ModelPart]' = {}
        for k, v in input_dict.items():
            _, intersected_model_part = ModelPartUtilities.IntersectModelParts(main_model_part, [v], add_neghbours)
            result[k] = intersected_model_part

        return result