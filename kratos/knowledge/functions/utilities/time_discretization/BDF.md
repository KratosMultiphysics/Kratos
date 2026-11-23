---
type: C++ Function
title: BDF
resource: utilities/time_discretization.h#L41-L44
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/time_discretization/BDF/SetAuxBDFPointer
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_time_discretization/TestTimeDiscretization/test_generic_BDF
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`BDF(const unsigned int TimeOrder) : mTimeOrder(TimeOrder)`

# Calls

- [SetAuxBDFPointer](../../../functions/utilities/time_discretization/BDF/SetAuxBDFPointer.md)

# Called by

- [test_generic_BDF](../../../functions/tests/test_time_discretization/TestTimeDiscretization/test_generic_BDF.md)