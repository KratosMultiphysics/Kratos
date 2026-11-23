---
type: C++ Function
title: ConvertToCsrMatrix
resource: utilities/amgcl_csr_conversion_utilities.h#L62-L82
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/SetIsOwnerOfData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/SetRowSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/SetColSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/AssignIndex1Data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/AssignIndex2Data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/AssignValueData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename CsrMatrix<TDataType, TIndexType>::UniquePointer ConvertToCsrMatrix( typename amgcl::backend::builtin<TDataType>::matrix& rA )`

# Calls

- [SetIsOwnerOfData](../../../functions/containers/csr_matrix/CsrMatrix/SetIsOwnerOfData.md)
- [SetRowSize](../../../functions/containers/csr_matrix/CsrMatrix/SetRowSize.md)
- [SetColSize](../../../functions/containers/csr_matrix/CsrMatrix/SetColSize.md)
- [AssignIndex1Data](../../../functions/containers/csr_matrix/CsrMatrix/AssignIndex1Data.md)
- [AssignIndex2Data](../../../functions/containers/csr_matrix/CsrMatrix/AssignIndex2Data.md)
- [AssignValueData](../../../functions/containers/csr_matrix/CsrMatrix/AssignValueData.md)