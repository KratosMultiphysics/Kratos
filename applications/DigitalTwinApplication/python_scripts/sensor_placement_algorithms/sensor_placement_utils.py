from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV

def AddSensorViewMasks(parameters: Kratos.Parameters, list_of_sensor_views: 'list[SensorViewUnionType]') -> None:
    masking_type = parameters["masking_type"].GetString()
    if masking_type == "automatic":
        for sensor_view in list_of_sensor_views:
            sensor_view.AddAuxiliaryExpression("mask", KratosDT.MaskUtils.GetMask(sensor_view.GetAuxiliaryExpression("filtered")))
    else:
        raise RuntimeError(f"Unsupported masking_type={masking_type} requested.")

def UpdateClusters(cluster_details: 'dict[int, tuple[list[int], ExpressionUnionType]]', list_of_sensor_views: 'list[SensorViewUnionType]') -> ExpressionUnionType:
    list_of_masks = [sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_sensor_views]
    cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks(list_of_masks)
    new_cluster_sensor_ids_list: 'list[list[int]]' = [[list_of_sensor_views[i].GetSensor().GetValue(KratosDT.SENSOR_ID) for i in indices_list] for indices_list, _ in cluster_data]

    list_of_cluster_ids_to_delete: 'list[int]' = list(cluster_details.keys())

    # first add ll new clusters
    for i, current_cluster_sensor_ids in enumerate(new_cluster_sensor_ids_list):
        old_cluster_id = -1
        for cluster_id, (cluster_sensor_ids, _) in cluster_details.items():
            if cluster_sensor_ids == current_cluster_sensor_ids:
                old_cluster_id = cluster_id
                del list_of_cluster_ids_to_delete[list_of_cluster_ids_to_delete.index(cluster_id)]
                break

        if old_cluster_id == -1:
            cluster_details[max(cluster_details.keys(), default=0) + 1] = (current_cluster_sensor_ids, cluster_data[i][1])
        else:
            cluster_details[old_cluster_id][1].SetExpression(cluster_data[i][1].GetExpression())

    # now remove clusters which are not anymore valid
    for k in list_of_cluster_ids_to_delete:
        del cluster_details[k]

    clustering_exp = list_of_masks[0] * 0.0
    for cluster_id, (_, cluster_mask_exp) in cluster_details.items():
        clustering_exp += cluster_mask_exp * cluster_id
    return clustering_exp

def OutputSensorViewExpressions(output_path: Path, vtu_output: Kratos.VtuOutput, sensor_view: SensorViewUnionType, additional_expressions: 'list[tuple[str, ExpressionUnionType]]' = []) -> None:
    vtu_output.ClearCellContainerExpressions()
    vtu_output.ClearNodalContainerExpressions()

    # add the main expression
    vtu_output.AddContainerExpression(sensor_view.GetExpressionName(), sensor_view.GetContainerExpression())

    # now add the auxiliary expression
    for auxiliary_name in sensor_view.GetAuxiliarySuffixes():
        vtu_output.AddContainerExpression(f"{sensor_view.GetExpressionName()}_{auxiliary_name}", sensor_view.GetAuxiliaryExpression(auxiliary_name))

    # now add the additional expressions
    for exp_name, exp in additional_expressions:
        vtu_output.AddContainerExpression(exp_name, exp)

    vtu_output.PrintOutput(str(output_path))

def PrintSensorViewsListToCSV(output_file_name: Path, list_of_sensor_views: 'list[SensorViewUnionType]', list_of_sensor_properties: 'list[str]') -> None:
    PrintSensorListToCSV(output_file_name, [sensor_view.GetSensor() for sensor_view in list_of_sensor_views], list_of_sensor_properties)

