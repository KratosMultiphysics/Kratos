---
type: C++ Method
title: WriteGeometryBlockMPI
resource: input_output/stl_io.cpp#L267-L317
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/stl_io/StlIO/IsValidGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteFacet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/stl_io/StlIO/WriteModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::WriteGeometryBlockMPI( const GeometryContainerType& rThisGeometries, const DataCommunicator& rDataCommunicator )`

# Calls

- [IsValidGeometry](../../../../functions/input_output/stl_io/StlIO/IsValidGeometry.md)
- [WriteFacet](../../../../functions/input_output/stl_io/StlIO/WriteFacet.md)
- [Recv](../../../../functions/includes/data_communicator/Recv.md)
- [Send](../../../../functions/includes/data_communicator/Send.md)

# Called by

- [WriteModelPart](../../../../functions/input_output/stl_io/StlIO/WriteModelPart.md)