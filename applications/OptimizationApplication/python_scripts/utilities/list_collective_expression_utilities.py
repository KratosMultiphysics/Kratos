import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA


def CollectiveListCollectiveProduct(collective_list: 'list[KratosOA.CollectiveExpression]', other_collective: KratosOA.CollectiveExpression) -> Kratos.Vector:
    result = Kratos.Vector(len(collective_list))

    for i, collective_list_item in enumerate(collective_list):
        result[i] = KratosOA.ExpressionUtils.InnerProduct(collective_list_item, other_collective)
    return result

def CollectiveListVectorProduct(collective_list: 'list[KratosOA.CollectiveExpression]', vector: Kratos.Vector) -> KratosOA.CollectiveExpression:
    if len(collective_list) != vector.Size():
        raise RuntimeError(f"Collective list size and vector size mismatch. [ Collective list size = {len(collective_list)}, vector size = {vector.Size()}]")
    if len(collective_list) == 0:
        raise RuntimeError("Collective lists cannot be empty.")

    result = collective_list[0] * 0.0
    for i, collective_list_item in enumerate(collective_list):
        result += collective_list_item * vector[i]

    return result