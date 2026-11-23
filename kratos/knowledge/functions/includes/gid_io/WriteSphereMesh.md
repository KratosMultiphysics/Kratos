---
type: C++ Function
title: WriteSphereMesh
resource: includes/gid_io.h#L1166-L1204
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_io_to_python/AddIOToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteSphereMesh(const MeshType& rThisMesh)`

# Calls

- [X0](../../../functions/includes/node/X0.md)
- [Y0](../../../functions/includes/node/Y0.md)
- [Z0](../../../functions/includes/node/Z0.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [AddIOToPython](../../../functions/python/add_io_to_python/AddIOToPython.md)