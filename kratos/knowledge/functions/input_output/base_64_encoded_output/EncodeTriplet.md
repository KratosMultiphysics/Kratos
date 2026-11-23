---
type: C++ Function
title: EncodeTriplet
resource: input_output/base_64_encoded_output.h#L291-L306
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/base_64_encoded_output/WriteData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void EncodeTriplet( std::ostream& rOutput, const std::array<char, 3>& rBytes, const IndexType Padding )`

# Calls

- [fill](../../../functions/includes/amatrix_interface/Matrix/fill.md)

# Called by

- [WriteData](../../../functions/input_output/base_64_encoded_output/WriteData.md)