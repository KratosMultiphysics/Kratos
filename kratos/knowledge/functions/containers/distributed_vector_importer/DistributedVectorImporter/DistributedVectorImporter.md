---
type: C++ Method
title: DistributedVectorImporter
resource: containers/distributed_vector_importer.h#L129-L137
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit DistributedVectorImporter(DistributedVectorImporter const& rOther) : mrComm(rOther.mrComm), mpNumbering(make_unique<DistributedNumbering<IndexType>>(*rOther.mpNumbering)), mLocalIByColor(rOther.mLocalIByColor), mLocallyOwnedIds(rOther.mLocallyOwnedIds), mIdOfLocallyOwnedTerms(rOther.mIdOfLocallyOwnedTerms), mVectorCommColors(rOther.mVectorCommColors)`