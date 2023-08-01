from typing import Any
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class ModelPartUtilities:
    @staticmethod
    def __GenerateUniqueModelPartName(prefix: str, list_of_names: 'list[str]', add_neighbours: bool) -> str:
        sorted_names = sorted(list_of_names)
        post_fix = "_IN" if add_neighbours else "_EN"
        return (f"{prefix}_" + "_".join(sorted_names) + post_fix).replace(".", "#")

    @staticmethod
    def UnionModelParts(main_model_part: Kratos.ModelPart, union_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> 'tuple[bool, Kratos.ModelPart]':
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName("Union", [model_part.FullName() for model_part in union_model_parts], add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return False, main_model_part.GetSubModelPart(unique_name)
        else:
            return True, Kratos.ModelPartOperationUtilities.Union(unique_name, main_model_part, union_model_parts, add_neighbours)

    @staticmethod
    def IntersectModelParts(main_model_part: Kratos.ModelPart, intersecting_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> 'tuple[bool, Kratos.ModelPart]':
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName("Intersect", [model_part.FullName() for model_part in intersecting_model_parts], add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return False, main_model_part.GetSubModelPart(unique_name)
        else:
            return True, Kratos.ModelPartOperationUtilities.Intersect(unique_name, main_model_part, intersecting_model_parts, add_neighbours)

    @staticmethod
    def SubstractModelParts(main_model_part: Kratos.ModelPart, substracting_model_parts: 'list[Kratos.ModelPart]', add_neighbours: bool) -> 'tuple[bool, Kratos.ModelPart]':
        unique_name = ModelPartUtilities.__GenerateUniqueModelPartName("Substract", [model_part.FullName() for model_part in substracting_model_parts], add_neighbours)
        if main_model_part.HasSubModelPart(unique_name):
            return False, main_model_part.GetSubModelPart(unique_name)
        else:
            return True, Kratos.ModelPartOperationUtilities.Substract(unique_name, main_model_part, substracting_model_parts, add_neighbours)

    @staticmethod
    def GetMergedMap(main_model_part: Kratos.ModelPart, input_dict: 'dict[Any, KratosOA.ContainerExpression.CollectiveExpressions]', add_neghbours: bool) -> 'dict[Any, Kratos.ModelPart]':
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