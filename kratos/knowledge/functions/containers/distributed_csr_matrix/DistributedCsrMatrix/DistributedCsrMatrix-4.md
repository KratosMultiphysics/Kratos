---
type: C++ Method
title: DistributedCsrMatrix
resource: containers/distributed_csr_matrix.h#L213-L229
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit DistributedCsrMatrix(const DistributedCsrMatrix& rOtherMatrix) : mpComm(rOtherMatrix.mpComm), mpRowNumbering(Kratos::make_unique< DistributedNumbering<TIndexType> >( rOtherMatrix.GetRowNumbering())), mpColNumbering(Kratos::make_unique< DistributedNumbering<TIndexType> >( rOtherMatrix.GetColNumbering())), mpDiagonalBlock(Kratos::make_unique<CsrMatrix<TDataType,TIndexType>>(rOtherMatrix.GetDiagonalBlock())), mpOffDiagonalBlock(Kratos::make_unique<CsrMatrix<TDataType,TIndexType>>(rOtherMatrix.GetOffDiagonalBlock())), mNonLocalData(rOtherMatrix.mNonLocalData), mOffDiagonalLocalIds(rOtherMatrix.mOffDiagonalLocalIds), mOffDiagonalGlobalIds(rOtherMatrix.mOffDiagonalGlobalIds), mfem_assemble_colors(rOtherMatrix.mfem_assemble_colors), mRecvCachedIJ(rOtherMatrix.mRecvCachedIJ), mSendCachedIJ(rOtherMatrix.mSendCachedIJ), mpVectorImporter(Kratos::make_unique<DistributedVectorImporter<TDataType,TIndexType>>(*rOtherMatrix.mpVectorImporter))`