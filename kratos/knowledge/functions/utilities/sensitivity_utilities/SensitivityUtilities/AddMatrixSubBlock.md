---
type: C++ Method
title: AddMatrixSubBlock
resource: utilities/sensitivity_utilities.cpp#L31-L63
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

`void SensitivityUtilities::AddMatrixSubBlock( Matrix& rOutput, const Matrix& rInput, const int RowOffset, const int ColOffset)`

# Called by

- [AssignEntityDerivativesToNodes](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes.md)