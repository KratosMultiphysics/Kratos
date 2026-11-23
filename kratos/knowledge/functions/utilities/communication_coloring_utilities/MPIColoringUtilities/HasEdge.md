---
type: C++ Method
title: HasEdge
resource: utilities/communication_coloring_utilities.cpp#L133-L142
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MPIColoringUtilities::HasEdge(std::map<int, std::map<int, int> >& rGraph, int i, int j)`

# Called by

- [ComputeCommunicationScheduling](../../../../functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling.md)