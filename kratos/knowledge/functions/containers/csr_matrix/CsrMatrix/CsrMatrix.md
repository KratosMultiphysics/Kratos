---
type: C++ Method
title: CsrMatrix
resource: containers/csr_matrix.h#L82-L85
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/ComputeColSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/ResizeValueData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/ResizeIndex1Data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/ResizeIndex2Data
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CsrMatrix() //needs to be public, since one could use the low level API to construc the CSR matrix`

# Calls

- [ComputeColSize](../../../../functions/containers/csr_matrix/CsrMatrix/ComputeColSize.md)
- [ResizeValueData](../../../../functions/containers/csr_matrix/CsrMatrix/ResizeValueData.md)
- [ResizeIndex1Data](../../../../functions/containers/csr_matrix/CsrMatrix/ResizeIndex1Data.md)
- [ResizeIndex2Data](../../../../functions/containers/csr_matrix/CsrMatrix/ResizeIndex2Data.md)