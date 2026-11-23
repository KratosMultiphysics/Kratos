---
type: C++ Method
title: CalculateEdgeDataValues
resource: containers/edge_based_data_structure.h#L449-L586
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetEdgeShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/SetLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/AddOffDiagonalValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetEdgeFaceShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateUnitNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/AddConvectiveBoundaryValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEdgeDataValues( const ModelPart& rModelPart, const SparseGraphType& rEdgesSparseGraph)`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)
- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetEdgeShapeFunctionsValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetEdgeShapeFunctionsValues.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetLength](../../../../functions/containers/edge_based_data_structure/EdgeData/SetLength.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [AddOffDiagonalValues](../../../../functions/containers/edge_based_data_structure/EdgeData/AddOffDiagonalValues.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [GetEdgeFaceShapeFunctionsValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetEdgeFaceShapeFunctionsValues.md)
- [CalculateUnitNormal](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateUnitNormal.md)
- [AddConvectiveBoundaryValue](../../../../functions/containers/edge_based_data_structure/EdgeData/AddConvectiveBoundaryValue.md)

# Called by

- [CalculateEdgeDataStructure](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataStructure.md)