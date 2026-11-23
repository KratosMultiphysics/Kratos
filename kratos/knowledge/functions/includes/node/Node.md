---
type: C++ Function
title: Node
resource: includes/node.h#L106-L116
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/CreateSolutionStepData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetNodalData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Node() : BaseType() , Flags() , mNodalData(0) , mDofs() , mData() , mInitialPosition() , mNodeLock()`

# Calls

- [CreateSolutionStepData](../../../functions/includes/node/CreateSolutionStepData.md)
- [SetNodalData](../../../functions/includes/dof/Dof/SetNodalData.md)