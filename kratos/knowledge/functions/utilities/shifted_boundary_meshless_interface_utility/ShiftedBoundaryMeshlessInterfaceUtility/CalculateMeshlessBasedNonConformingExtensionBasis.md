---
type: C++ Method
title: CalculateMeshlessBasedNonConformingExtensionBasis
resource: utilities/shifted_boundary_meshless_interface_utility.cpp#L679-L768
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags
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
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetMLSShapeFunctionsAndGradientsFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSplitElementSupportCloud
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
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateKernelRadius
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ShiftedBoundaryMeshlessInterfaceUtility::CalculateMeshlessBasedNonConformingExtensionBasis()`

# Calls

- [SetInterfaceFlags](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetStandardModifiedShapeFunctionsFactory](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/GetStandardModifiedShapeFunctionsFactory.md)
- [GetMLSShapeFunctionsAndGradientsFunction](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetMLSShapeFunctionsAndGradientsFunction.md)
- [GetElementSizeFunction](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/GetElementSizeFunction.md)
- [SetSplitElementSupportCloud](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSplitElementSupportCloud.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [SetNodalDistancesVector](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/SetNodalDistancesVector.md)
- [ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues.md)
- [ComputePositiveSideInterfaceAreaNormals](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputePositiveSideInterfaceAreaNormals.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [CalculateKernelRadius](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateKernelRadius.md)

# Called by

- [CalculateExtensionOperator](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateExtensionOperator.md)