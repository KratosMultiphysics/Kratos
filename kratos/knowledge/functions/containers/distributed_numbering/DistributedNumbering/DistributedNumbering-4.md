---
type: C++ Method
title: DistributedNumbering
resource: containers/distributed_numbering.h#L101-L119
generated:
  by: okf-rs/0.7.0
---

# Signature

`DistributedNumbering( const DataCommunicator* pComm, const IndexType TotalSize, const MpiIndexType Nranks //note that we could obtain this by the rComm, but we need to distinguish this constructor from the previous ) : mpComm(pComm)`