---
type: C++ Module
title: calculate_signed_distance_to_3d_skin_process
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1-L2848
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/ctime
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-deprecated-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spatial-containers-octree-binary-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-spatial-containers-configure-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-math-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometry-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-3d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-body-normal-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-flags-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-binbased-fast-point-locator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-binbased-nodes-in-element-locator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-calculate-distance-to-skin-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/omp-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistanceSpatialContainersConfigure](../../classes/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure.md)
- [CellNodeData](../../classes/processes/calculate_signed_distance_to_3d_skin_process/CellNodeData.md)
- [DistanceSpatialContainersConfigure](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/DistanceSpatialContainersConfigure.md)
- [AllocateData](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/AllocateData.md)
- [CopyData](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/CopyData.md)
- [DeleteData](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/DeleteData.md)
- [CalculateBoundingBox](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/CalculateBoundingBox.md)
- [GetBoundingBox](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/GetBoundingBox.md)
- [Intersection](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/Intersection.md)
- [IntersectionBox](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/IntersectionBox.md)
- [IsIntersected](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/IsIntersected.md)
- [Info](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/Info.md)
- [PrintInfo](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/PrintInfo.md)
- [PrintData](../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/PrintData.md)
- [CalculateSignedDistanceTo3DSkinProcess](../../classes/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess.md)
- [CalculateSignedDistanceTo3DSkinProcess](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateSignedDistanceTo3DSkinProcess.md)
- [Execute](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/Execute.md)
- [MappingPressureToStructure](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure.md)
- [ComputeDiscontinuousInterpolation](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation.md)
- [AveragePressureToNode](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AveragePressureToNode.md)
- [DistanceFluidStructure](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/DistanceFluidStructure.md)
- [Initialize](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/Initialize.md)
- [SetIndexTable](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/SetIndexTable.md)
- [CalcElementDistances](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcElementDistances.md)
- [IdentifyIntersectionNodes](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IdentifyIntersectionNodes.md)
- [StructuralElementNotYetConsidered](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/StructuralElementNotYetConsidered.md)
- [IsIntersectionNodeOnTetEdge](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsIntersectionNodeOnTetEdge.md)
- [IsNewIntersectionNode](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsNewIntersectionNode.md)
- [IsIntersectionOnCorner](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsIntersectionOnCorner.md)
- [CalculateNormal3D](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNormal3D.md)
- [CalcDistanceTo3DSkin](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcDistanceTo3DSkin.md)
- [FillIntNodesContainer](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/FillIntNodesContainer.md)
- [CalcSignedDistancesToOneIntNode](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToOneIntNode.md)
- [CalcSignedDistancesToTwoIntNodes](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToTwoIntNodes.md)
- [CalcSignedDistancesToThreeIntNodes](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToThreeIntNodes.md)
- [CalcSignedDistancesToMoreThanThreeIntNodes](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)
- [PointDistanceToPlane](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PointDistanceToPlane.md)
- [AssignMinimalNodalDistance](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AssignMinimalNodalDistance.md)
- [AvoidZeroDistances](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AvoidZeroDistances.md)
- [GenerateSkinModelPart](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart.md)
- [GenerateOctree](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateOctree.md)
- [GenerateNodes](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateNodes.md)
- [GenerateCellNode](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateCellNode.md)
- [SetNodeInNeighbours](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/SetNodeInNeighbours.md)
- [CalculateDistance2](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance2.md)
- [CalculateDistance](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance.md)
- [CalculateDistance](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance-2.md)
- [CalculateNotEmptyLeavesDistance](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [CalculateNodeDistance](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNodeDistance.md)
- [DistancePositionInSpace](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/DistancePositionInSpace.md)
- [GetIntersectionsAndNodes](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes.md)
- [GetIntersections](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersections.md)
- [GetCellIntersections](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetCellIntersections.md)
- [IntersectionTriangleSegment](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IntersectionTriangleSegment.md)
- [Info](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/Info.md)
- [PrintInfo](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintInfo.md)
- [PrintData](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintData.md)
- [PrintGiDMesh](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintGiDMesh.md)
- [PrintGiDResults](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/PrintGiDResults.md)
- [EigenVectors](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/EigenVectors.md)
- [CreatePartition](../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CreatePartition.md)

# Imports

- `<string>`
- `<iostream>`
- `<ctime>`
- `"includes/define.h"`
- `"processes/process.h"`
- `"includes/model_part.h"`
- `"includes/deprecated_variables.h"`
- `"spatial_containers/octree_binary.h"`
- `"utilities/spatial_containers_configure.h"`
- `"utilities/timer.h"`
- `"utilities/math_utils.h"`
- `"utilities/geometry_utilities.h"`
- `"geometries/triangle_3d_3.h"`
- `"geometries/quadrilateral_3d_4.h"`
- `"utilities/body_normal_calculation_utils.h"`
- `"includes/kratos_flags.h"`
- `"utilities/binbased_fast_point_locator.h"`
- `"utilities/binbased_nodes_in_element_locator.h"`
- `"processes/calculate_distance_to_skin_process.h"`
- `"omp.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)