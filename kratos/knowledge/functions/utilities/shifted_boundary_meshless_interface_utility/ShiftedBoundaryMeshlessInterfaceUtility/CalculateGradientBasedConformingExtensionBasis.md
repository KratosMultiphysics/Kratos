---
type: C++ Method
title: CalculateGradientBasedConformingExtensionBasis
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L212-L482
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSurrogateBoundaryNodalGradientWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/GetStandardModifiedShapeFunctionsFactory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/SetNodalDistancesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ShiftedBoundaryMeshlessInterfaceUtility::CalculateGradientBasedConformingExtensionBasis()`

# Calls

- [SetInterfaceFlags](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags.md)
- [SetSurrogateBoundaryNodalGradientWeights](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSurrogateBoundaryNodalGradientWeights.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetStandardModifiedShapeFunctionsFactory](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/GetStandardModifiedShapeFunctionsFactory.md)
- [GetElementSizeFunction](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [SetNodalDistancesVector](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/SetNodalDistancesVector.md)
- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputePositiveSideInterfaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalculateExtensionOperator](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator.md)