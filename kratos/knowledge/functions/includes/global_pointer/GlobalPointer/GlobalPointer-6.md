---
type: C++ Method
title: GlobalPointer
resource: includes/global_pointer.h#L128-L139
generated:
  by: okf-rs/0.7.0
---

# Signature

`GlobalPointer( Kratos::weak_ptr<TDataType> DataPointer, int Rank = 0 ) : mDataPointer(DataPointer.lock().get()) #ifdef KRATOS_USING_MPI , mRank(Rank) #endif`