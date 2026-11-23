---
type: C++ Function
title: WriteNodalResults
resource: includes/gid_io.h#L642-L660
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

`void WriteNodalResults( Variable<bool> const& rVariable, const NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber )`

# Calls

- [Start](../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)