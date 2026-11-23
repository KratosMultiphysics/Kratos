---
type: C++ Function
title: Node
resource: includes/node.h#L272-L282
generated:
  by: okf-rs/0.7.0
---

# Signature

`Node(IndexType NewId, const double NewX, const double NewY, const double NewZ, VariablesList::Pointer pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1) : BaseType(NewX, NewY, NewZ) , Flags() , mNodalData(NewId, pVariablesList,ThisData,NewQueueSize) , mDofs() , mData() , mInitialPosition(NewX, NewY, NewZ) , mNodeLock()`