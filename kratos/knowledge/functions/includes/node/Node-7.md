---
type: C++ Function
title: Node
resource: includes/node.h#L207-L223
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit Node(const Node& rOtherNode) : BaseType(rOtherNode) , Flags(rOtherNode) , mNodalData(rOtherNode.Id()) , mData(rOtherNode.mData) , mInitialPosition(rOtherNode.mInitialPosition) , mNodeLock()`