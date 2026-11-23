---
type: C++ Method
title: AssignIndex1Data
resource: containers/csr_matrix.h#L331-L340
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/IsOwnerOfData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/csr_matrix/CsrMatrix/Clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/amgcl_csr_conversion_utilities/ConvertToCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignIndex1Data(TIndexType* pExternalData, TIndexType DataSize)`

# Calls

- [IsOwnerOfData](../../../../functions/containers/csr_matrix/CsrMatrix/IsOwnerOfData.md)

# Called by

- [Clear](../../../../functions/containers/csr_matrix/CsrMatrix/Clear.md)
- [ConvertToCsrMatrix](../../../../functions/utilities/amgcl_csr_conversion_utilities/ConvertToCsrMatrix.md)