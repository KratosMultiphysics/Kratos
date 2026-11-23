---
type: C++ Module
title: calculate_signed_distance_to_3d_condition_skin_process
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1-L2176
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
  - target: external/utilities-body-normal-calculation-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [DistanceSpatialContainersConditionConfigure](../../classes/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure.md)
- [CellNodeData](../../classes/processes/calculate_signed_distance_to_3d_condition_skin_process/CellNodeData.md)
- [DistanceSpatialContainersConditionConfigure](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/DistanceSpatialContainersConditionConfigure.md)
- [AllocateData](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/AllocateData.md)
- [CopyData](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/CopyData.md)
- [DeleteData](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/DeleteData.md)
- [CalculateBoundingBox](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/CalculateBoundingBox.md)
- [GetBoundingBox](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/GetBoundingBox.md)
- [Intersection](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/Intersection.md)
- [IntersectionBox](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/IntersectionBox.md)
- [IsIntersected](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/IsIntersected.md)
- [Info](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/Info.md)
- [PrintInfo](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/PrintInfo.md)
- [PrintData](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/PrintData.md)
- [CalculateSignedDistanceTo3DConditionSkinProcess](../../classes/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess.md)
- [CalculateSignedDistanceTo3DConditionSkinProcess](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateSignedDistanceTo3DConditionSkinProcess.md)
- [Execute](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/Execute.md)
- [DistanceFluidStructure](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/DistanceFluidStructure.md)
- [InitializeDistances](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/InitializeDistances.md)
- [SetIndexTable](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/SetIndexTable.md)
- [CalcNodalDistancesOfTetNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistancesOfTetNodes.md)
- [IdentifyIntersectionNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IdentifyIntersectionNodes.md)
- [StructuralElementNotYetConsidered](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/StructuralElementNotYetConsidered.md)
- [IsIntersectionNodeOnTetEdge](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsIntersectionNodeOnTetEdge.md)
- [IsNewIntersectionNode](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsNewIntersectionNode.md)
- [IsIntersectionOnCorner](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsIntersectionOnCorner.md)
- [CalculateNormal3D](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNormal3D.md)
- [CalcNodalDistanceTo3DSkin](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin.md)
- [ComputeApproximationNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/ComputeApproximationNodes.md)
- [CalcSignedDistancesToOneIntNode](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToOneIntNode.md)
- [CalcSignedDistancesToTwoIntNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToTwoIntNodes.md)
- [CalcSignedDistancesToThreeIntNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToThreeIntNodes.md)
- [CalcSignedDistancesToFourIntNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToFourIntNodes.md)
- [FindIndexNodesOfTriangle2](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/FindIndexNodesOfTriangle2.md)
- [CalcSignedDistancesToApproxTriangle](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle.md)
- [AssignDistancesToElements](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/AssignDistancesToElements.md)
- [GenerateOctree](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateOctree.md)
- [GenerateNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateNodes.md)
- [GenerateCellNode](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateCellNode.md)
- [SetNodeInNeighbours](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/SetNodeInNeighbours.md)
- [CalculateDistance2](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance2.md)
- [CalculateDistance](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance.md)
- [CalculateDistance](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance-2.md)
- [CalculateNotEmptyLeavesDistance](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [CalculateNodeDistance](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNodeDistance.md)
- [DistancePositionInSpace](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/DistancePositionInSpace.md)
- [GetIntersectionsAndNodes](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes.md)
- [GetIntersections](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersections.md)
- [GetCellIntersections](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetCellIntersections.md)
- [IntersectionTriangleSegment](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IntersectionTriangleSegment.md)
- [Info](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/Info.md)
- [PrintInfo](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintInfo.md)
- [PrintData](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintData.md)
- [PrintGiDMesh](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintGiDMesh.md)
- [PrintGiDResults](../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/PrintGiDResults.md)

# Imports

- `<string>`
- `<iostream>`
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
- `"utilities/body_normal_calculation_utils.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)