---
type: C++ Function
title: AddConnectivityData
resource: input_output/vtu_output.cpp#L307-L337
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetGeometryTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetOffsets
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetWritingNDData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddConnectivityData( XmlElementsArray& rCellElement, std::vector<IndexType>& rWritingIndices, const ModelPart::NodesContainerType& rNodes, VtuOutput::CellContainerPointerType pCells, TXmlDataElementWrapper& rXmlDataElementWrapper, const IndexType EchoLevel)`

# Calls

- [visit](../../../functions/python/variable_indexing_python/visit.md)
- [GetGeometryTypes](../../../functions/input_output/vtu_output/GetGeometryTypes.md)
- [GetOffsets](../../../functions/input_output/vtu_output/GetOffsets.md)
- [GetConnectivities](../../../functions/input_output/vtu_output/GetConnectivities.md)
- [GetWritingNDData](../../../functions/input_output/vtu_output/GetWritingNDData.md)

# Called by

- [WriteUnstructuredGridData](../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)