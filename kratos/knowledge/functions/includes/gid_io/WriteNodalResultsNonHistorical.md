---
type: C++ Function
title: WriteNodalResultsNonHistorical
resource: includes/gid_io.h#L870-L885
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

`void WriteNodalResultsNonHistorical( Variable<bool> const& rVariable, const NodesContainerType& rNodes, const double SolutionTag )`

# Calls

- [Start](../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)