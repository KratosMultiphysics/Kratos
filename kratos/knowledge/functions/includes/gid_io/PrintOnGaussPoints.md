---
type: C++ Function
title: PrintOnGaussPoints
resource: includes/gid_io.h#L1363-L1379
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void PrintOnGaussPoints( const Variable<bool>& rVariable, const ModelPart& rModelPart, const double SolutionTag, const int ValueIndex = 0 )`

# Calls

- [Start](../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)