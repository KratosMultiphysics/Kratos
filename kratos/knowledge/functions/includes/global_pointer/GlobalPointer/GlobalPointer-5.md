---
type: C++ Method
title: GlobalPointer
resource: includes/global_pointer.h#L111-L122
generated:
  by: okf-rs/0.7.0
---

# Signature

`GlobalPointer( Kratos::intrusive_ptr<TDataType>& DataPointer, int Rank = 0 ) : mDataPointer(DataPointer.get()) #ifdef KRATOS_USING_MPI , mRank(Rank) #endif`