---
type: C++ Method
title: SetAuxBDFPointer
resource: utilities/time_discretization.cpp#L29-L56
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/time_discretization/BDF
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BDF::SetAuxBDFPointer( const std::size_t TimeOrder, Kratos::unique_ptr<BDF> &rpAuxBDF)`

# Called by

- [BDF](../../../../functions/utilities/time_discretization/BDF.md)