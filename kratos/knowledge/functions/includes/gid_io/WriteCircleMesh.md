---
type: C++ Function
title: WriteCircleMesh
resource: includes/gid_io.h#L1206-L1246
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
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_io_to_python/AddIOToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteCircleMesh(const MeshType& rThisMesh)`

# Calls

- [X0](../../../functions/includes/node/X0.md)
- [Y0](../../../functions/includes/node/Y0.md)
- [Z0](../../../functions/includes/node/Z0.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [AddIOToPython](../../../functions/python/add_io_to_python/AddIOToPython.md)