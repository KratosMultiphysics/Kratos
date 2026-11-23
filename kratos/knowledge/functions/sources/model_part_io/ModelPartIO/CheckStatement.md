---
type: C++ Method
title: CheckStatement
resource: sources/model_part_io.cpp#L6552-L6555
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::CheckStatement(std::string const& rStatement, std::string const& rGivenWord) const`

# Called by

- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)