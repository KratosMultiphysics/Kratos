---
type: C++ Method
title: ResizeIndex2Data
resource: containers/csr_matrix.h#L373-L380
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/IsOwnerOfData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/csr_matrix/CsrMatrix/CsrMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/reserve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ResizeIndex2Data(TIndexType DataSize)`

# Calls

- [IsOwnerOfData](../../../../functions/containers/csr_matrix/CsrMatrix/IsOwnerOfData.md)

# Called by

- [CsrMatrix](../../../../functions/containers/csr_matrix/CsrMatrix/CsrMatrix.md)
- [reserve](../../../../functions/containers/csr_matrix/CsrMatrix/reserve.md)