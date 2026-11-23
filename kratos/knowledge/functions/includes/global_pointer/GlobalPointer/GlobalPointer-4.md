---
type: C++ Method
title: GlobalPointer
resource: includes/global_pointer.h#L94-L105
generated:
  by: okf-rs/0.7.0
---

# Signature

`GlobalPointer( Kratos::shared_ptr<TDataType> DataPointer, int Rank = 0 ) : mDataPointer(DataPointer.get()) #ifdef KRATOS_USING_MPI , mRank(Rank) #endif`