---
type: C++ Method
title: GetMatrixSubBlock
resource: utilities/sensitivity_utilities.cpp#L65-L102
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SensitivityUtilities::GetMatrixSubBlock( Matrix& rOutput, const Matrix& rInput, const int RowOffset, const int Rows, const int ColOffset, const int Cols)`

# Called by

- [AssignEntityDerivativesToNodes](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes.md)