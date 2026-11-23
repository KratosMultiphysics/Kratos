---
type: C++ Method
title: NodalData
resource: containers/nodal_data.cpp#L40-L45
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`NodalData::NodalData(IndexType TheId, VariablesList::Pointer pVariablesList, BlockType const * ThisData, SizeType NewQueueSize) : mId(TheId) , mSolutionStepsNodalData(pVariablesList,ThisData,NewQueueSize)`