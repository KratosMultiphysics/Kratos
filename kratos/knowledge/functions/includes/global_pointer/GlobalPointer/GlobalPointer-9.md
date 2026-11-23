---
type: C++ Method
title: GlobalPointer
resource: includes/global_pointer.h#L163-L169
generated:
  by: okf-rs/0.7.0
---

# Signature

`GlobalPointer(const GlobalPointer&& rOther) : mDataPointer(std::move(rOther.mDataPointer)) #ifdef KRATOS_USING_MPI , mRank(std::move(rOther.mRank)) #endif`