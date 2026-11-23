---
type: C++ Method
title: ExecuteInitialize
resource: processes/flux_corrected_transport_convection_process.h#L124-L149
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExecuteInitialize() override`

# Calls

- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)
- [CalculateEdgeDataStructure](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataStructure.md)