---
type: C++ Method
title: ResultsProxy
resource: utilities/pointer_communicator.h#L81-L88
generated:
  by: okf-rs/0.7.0
---

# Signature

`ResultsProxy( int current_rank, GlobalPointersUnorderedMap< TPointerDataType, TSendType > NonLocalData, TFunctorType UserFunctor, GlobalPointerCommunicator<TPointerDataType>* pPointerComm ): mCurrentRank(current_rank), mNonLocalData(NonLocalData), mUserFunctor(UserFunctor), mpPointerComm(pPointerComm)`