---
type: C++ Function
title: ComputeCommunicationPlan
resource: utilities/pointer_map_communicator.h#L478-L489
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<int> ComputeCommunicationPlan(const std::unordered_map<int, TArgs...>& rNonLocalGlobalPointerMap)`

# Calls

- [ComputeCommunicationScheduling](../../../functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling.md)