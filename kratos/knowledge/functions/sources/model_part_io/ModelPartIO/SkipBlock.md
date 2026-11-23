---
type: C++ Method
title: SkipBlock
resource: sources/model_part_io.cpp#L1562-L1589
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckStatement
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodesNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadInitialValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadMeshBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::SkipBlock(std::string const& BlockName)`

# Calls

- [CheckStatement](../../../../functions/sources/model_part_io/ModelPartIO/CheckStatement.md)

# Called by

- [ReadNodes](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodes.md)
- [ReadNodesNumber](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodesNumber.md)
- [ReadProperties](../../../../functions/sources/model_part_io/ModelPartIO/ReadProperties.md)
- [ReadGeometries](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometries.md)
- [ReadGeometriesConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities.md)
- [ReadElements](../../../../functions/sources/model_part_io/ModelPartIO/ReadElements.md)
- [ReadElementsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities.md)
- [ReadConditions](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditions.md)
- [ReadConditionsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities.md)
- [ReadConstraints](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraints.md)
- [ReadConstraintsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities.md)
- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)
- [ReadNodalGraph](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)
- [ReadNodalGraphFromEntitiesList](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList.md)
- [ReadNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)
- [ReadCommunicatorDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock.md)
- [ReadMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadMeshBlock.md)
- [ReadSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock.md)
- [DivideMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshBlock.md)
- [DivideSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock.md)