---
type: C++ Function
title: WriteLocalAxesOnNodes
resource: includes/gid_io.h#L816-L838
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
  called_by:
  - target: functions/python/add_io_to_python/AddIOToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void WriteLocalAxesOnNodes( Variable<array_1d<double, 3> > const& rVariable, const NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber )`

# Calls

- [Start](../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [AddIOToPython](../../../functions/python/add_io_to_python/AddIOToPython.md)