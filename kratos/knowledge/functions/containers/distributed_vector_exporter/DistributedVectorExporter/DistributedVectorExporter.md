---
type: C++ Method
title: DistributedVectorExporter
resource: containers/distributed_vector_exporter.h#L126-L134
generated:
  by: okf-rs/0.7.0
---

# Signature

`DistributedVectorExporter(DistributedVectorExporter const& rOther) : mrComm(rOther.mrComm), mpNumbering(Kratos::make_unique<DistributedNumbering<TIndexType>>(*rOther.mpNumbering)), mToRecvLocalId(rOther.mToRecvLocalId), mPositionWithinData(rOther.mPositionWithinData), mVectorCommColors(rOther.mVectorCommColors)`