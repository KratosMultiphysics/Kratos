---
type: C++ Function
title: IsInParallel
resource: utilities/openmp_utils.h#L122-L129
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/variables_list/AddDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Free
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/pointer_map_communicator/ApplyProxy/ApplyProxy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/pointer_map_communicator/if
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline int IsInParallel()`

# Called by

- [AddDof](../../../functions/containers/variables_list/AddDof.md)
- [Fix](../../../functions/includes/node/Fix.md)
- [Free](../../../functions/includes/node/Free.md)
- [ApplyProxy](../../../functions/utilities/pointer_map_communicator/ApplyProxy/ApplyProxy.md)
- [if](../../../functions/utilities/pointer_map_communicator/if.md)