//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_DISTRIBUTED_CSR_MATRIX_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_CSR_MATRIX_H_INCLUDED


// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/csr_matrix.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_vector_importer.h"
#include "containers/distributed_vector_exporter.h"
#include "includes/key_hash.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// This class implements "serial" CSR matrix, including capabilities for FEM assembly
template< class TDataType=double, class TIndexType=std::size_t>
class DistributedCsrMatrix final
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    using SizeType = IndexType;
    typedef int MpiIndexType;
    typedef CsrMatrix<TDataType,TIndexType> BlockMatrixType;
    typedef typename CsrMatrix<TDataType,TIndexType>::MatrixMapType MatrixMapType;

    /// Pointer definition of DistributedCsrMatrix
    KRATOS_CLASS_POINTER_DEFINITION(DistributedCsrMatrix);

    ///@}
    ///@name Life Cycle
    ///@{

    //default constructor. To be used for low level operations
    DistributedCsrMatrix()
    {}

    DistributedCsrMatrix(const DataCommunicator& rComm)
        : mpComm(&rComm)
    {}

    DistributedCsrMatrix(const DistributedSparseGraph<TIndexType>& rSparseGraph)
        :
        mpComm(&rSparseGraph.GetComm())
    {
        mpRowNumbering = Kratos::make_unique< DistributedNumbering<TIndexType> >( rSparseGraph.GetRowNumbering());

        //compute the columns size
        TIndexType max_col_index = rSparseGraph.ComputeMaxGlobalColumnIndex();
        TIndexType tot_col_size = max_col_index+1;

        //this ensures that diagonal blocks are square for square matrices
        if (tot_col_size == mpRowNumbering->Size()) {
            mpColNumbering = Kratos::make_unique<DistributedNumbering<TIndexType>>(*mpComm, mpRowNumbering->GetCpuBounds());
        } else {
            mpColNumbering = Kratos::make_unique<DistributedNumbering<TIndexType>>(*mpComm, tot_col_size, mpComm->Size());
        }

        mOffDiagonalLocalIds.clear(); //this is the map that allows to transform from global_ids to local Ids for entries in the non_diag block

        //count entries in diagonal and off_diag blocks
        const auto& local_graph = rSparseGraph.GetLocalGraph();
        const TIndexType nlocal_rows = local_graph.Size();

        TIndexType diag_nnz = 0, offdiag_nnz=0;
        for(const auto& entries : local_graph) {
            for(const auto& global_j : entries) {
                if(GetColNumbering().IsLocal(global_j)) {
                    diag_nnz += 1;
                } else {
                    offdiag_nnz += 1;
                    mOffDiagonalLocalIds[global_j] = 0;
                }
            }
        }

        //begin by computing local Ids for non diagonal block
        TIndexType counter = 0;
        for(auto& item : mOffDiagonalLocalIds)
            item.second = counter++;

        mOffDiagonalGlobalIds.resize(mOffDiagonalLocalIds.size(), false);
        counter = 0;
        for(auto& item : mOffDiagonalLocalIds) {
            TIndexType global_i = item.first;
            mOffDiagonalGlobalIds[counter++] = global_i;
        }

        //*******************************
        //construct diagonal block
        GetDiagonalBlock().reserve(nlocal_rows, diag_nnz);

        if(diag_nnz == 0) { //empty matrix
            for(unsigned int i=0; i<nlocal_rows+1; ++i)
                GetDiagonalBlock().index1_data()[i] = 0;
            GetDiagonalBlock().SetColSize(GetColNumbering().LocalSize());
        } else {
            counter = 0;
            for(const auto& entries : local_graph) {
                unsigned int k = 0;
                TIndexType row_begin = GetDiagonalBlock().index1_data()[counter];
                for(auto global_j : entries) {
                    if(GetColNumbering().IsLocal(global_j)) {
                        TIndexType local_j = GetColNumbering().LocalId(global_j);
                        GetDiagonalBlock().index2_data()[row_begin+k] = local_j;
                        GetDiagonalBlock().value_data()[row_begin+k] = 0.0;
                        k++;
                    }
                }
                GetDiagonalBlock().index1_data()[counter+1] = row_begin + k;
                counter++;
            }

            GetDiagonalBlock().SetColSize(GetColNumbering().LocalSize());
        }
#ifdef KRATOS_DEBUG
        GetDiagonalBlock().CheckColSize();
#endif
        //ensure columns are ordered
        for(TIndexType i = 0; i<nlocal_rows; ++i) {
            TIndexType row_begin = GetDiagonalBlock().index1_data()[i];
            TIndexType row_end = GetDiagonalBlock().index1_data()[i+1];

            if(row_end - row_begin > 0)
                std::sort(GetDiagonalBlock().index2_data().begin() + row_begin,GetDiagonalBlock().index2_data().begin()+row_end);
        }

        //*******************************
        //construct offdiagonal block

        //store off-diagonal block
        GetOffDiagonalBlock().reserve(nlocal_rows, offdiag_nnz);
        if(offdiag_nnz == 0) { //empty matrix
            for(unsigned int i=0; i<nlocal_rows+1; ++i)
                GetOffDiagonalBlock().index1_data()[i] = 0;
            GetOffDiagonalBlock().SetColSize(mOffDiagonalLocalIds.size());
        } else {
            counter = 0;
            for(const auto& entries : local_graph) {
                unsigned int k = 0;
                TIndexType row_begin = GetOffDiagonalBlock().index1_data()[counter];
                for(auto global_j : entries) {
                    if( ! GetColNumbering().IsLocal(global_j)) {
                        TIndexType local_j = GetOffDiagonalBlockLocalId(global_j);
                        GetOffDiagonalBlock().index2_data()[row_begin+k] = local_j;
                        GetOffDiagonalBlock().value_data()[row_begin+k] = 0.0;
                        k++;
                    }
                }
                GetOffDiagonalBlock().index1_data()[counter+1] = row_begin + k;
                counter++;
            }
            GetOffDiagonalBlock().SetColSize(mOffDiagonalLocalIds.size());
        }
#ifdef KRATOS_DEBUG
        GetOffDiagonalBlock().CheckColSize();
#endif

        //ensure columns are ordered
        for(TIndexType i = 0; i<nlocal_rows; ++i) {
            TIndexType row_begin = GetOffDiagonalBlock().index1_data()[i];
            TIndexType row_end = GetOffDiagonalBlock().index1_data()[i+1];
            if(row_end - row_begin > 0)
                std::sort(GetOffDiagonalBlock().index2_data().begin() + row_begin,GetOffDiagonalBlock().index2_data().begin()+row_end);
        }

        PrepareNonLocalCommunications(rSparseGraph);

        //mount importer for SpMV calculations
        auto pimporter = Kratos::make_unique<DistributedVectorImporter<TDataType,TIndexType>>(GetComm(),mOffDiagonalGlobalIds, GetColNumbering());
        mpVectorImporter.swap(pimporter);
    }

    explicit DistributedCsrMatrix(const DistributedCsrMatrix& rOtherMatrix)
        :
        mpComm(rOtherMatrix.mpComm),
        mpRowNumbering(Kratos::make_unique< DistributedNumbering<TIndexType> >( rOtherMatrix.GetRowNumbering())),
        mpColNumbering(Kratos::make_unique< DistributedNumbering<TIndexType> >( rOtherMatrix.GetColNumbering())),
        mpDiagonalBlock(Kratos::make_unique<CsrMatrix<TDataType,TIndexType>>(rOtherMatrix.GetDiagonalBlock())),
        mpOffDiagonalBlock(Kratos::make_unique<CsrMatrix<TDataType,TIndexType>>(rOtherMatrix.GetOffDiagonalBlock())),
        mNonLocalData(rOtherMatrix.mNonLocalData),
        mOffDiagonalLocalIds(rOtherMatrix.mOffDiagonalLocalIds),
        mOffDiagonalGlobalIds(rOtherMatrix.mOffDiagonalGlobalIds),
        mfem_assemble_colors(rOtherMatrix.mfem_assemble_colors),
        mRecvCachedIJ(rOtherMatrix.mRecvCachedIJ),
        mSendCachedIJ(rOtherMatrix.mSendCachedIJ),
        mpVectorImporter(Kratos::make_unique<DistributedVectorImporter<TDataType,TIndexType>>(*rOtherMatrix.mpVectorImporter))
    {
        ReconstructDirectAccessVectors();
    }

    //move constructor
    DistributedCsrMatrix(DistributedCsrMatrix<TDataType,TIndexType>&& rOtherMatrix)
        :mpComm(rOtherMatrix.mpComm)
    {
        mpRowNumbering.swap(rOtherMatrix.pGetRowNumbering());
        mpColNumbering.swap(rOtherMatrix.pGetColNumbering());
        mpDiagonalBlock = std::move(rOtherMatrix.mpDiagonalBlock);
        mpOffDiagonalBlock = std::move(rOtherMatrix.mpOffDiagonalBlock);
        mNonLocalData = std::move(rOtherMatrix.mNonLocalData);
        mSendCachedIJ = std::move(rOtherMatrix.mSendCachedIJ);
        mRecvCachedIJ = std::move(rOtherMatrix.mRecvCachedIJ);
        mOffDiagonalLocalIds = std::move(rOtherMatrix.mOffDiagonalLocalIds);
        mOffDiagonalGlobalIds = std::move(rOtherMatrix.mOffDiagonalGlobalIds);
        mfem_assemble_colors = std::move(rOtherMatrix.mfem_assemble_colors);
        mpVectorImporter = std::move(rOtherMatrix.mpVectorImporter);
    }

    //move assignement operator
    DistributedCsrMatrix& operator=(DistributedCsrMatrix&& rOtherMatrix)
    {
        mpComm = rOtherMatrix.mpComm,
        mpRowNumbering.swap(rOtherMatrix.pGetRowNumbering());
        mpColNumbering.swap(rOtherMatrix.pGetColNumbering());
        mpDiagonalBlock = std::move(rOtherMatrix.mpDiagonalBlock);
        mpOffDiagonalBlock = std::move(rOtherMatrix.mpOffDiagonalBlock);
        mNonLocalData = std::move(rOtherMatrix.mNonLocalData);
        mSendCachedIJ = std::move(rOtherMatrix.mSendCachedIJ);
        mRecvCachedIJ = std::move(rOtherMatrix.mRecvCachedIJ);
        mOffDiagonalLocalIds = std::move(rOtherMatrix.mOffDiagonalLocalIds);
        mOffDiagonalGlobalIds = std::move(rOtherMatrix.mOffDiagonalGlobalIds);
        mfem_assemble_colors = std::move(rOtherMatrix.mfem_assemble_colors);
        mpVectorImporter = std::move(rOtherMatrix.mpVectorImporter);
    }
    /// Destructor.
    ~DistributedCsrMatrix() {}

    /// Assignment operator. TODO: decide if we do want to allow it
    DistributedCsrMatrix& operator=(DistributedCsrMatrix const& rOther)=delete;
    // {
    //     this->AddEntries(rOther.GetGraph());
    //     return *this;
    // }

    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
    }

    inline const DistributedNumbering<TIndexType>& GetRowNumbering() const
    {
        return *mpRowNumbering;
    }

    inline const DistributedNumbering<TIndexType>& GetColNumbering() const
    {
        return *mpColNumbering;
    }


    inline typename DistributedNumbering<TIndexType>::UniquePointer& pGetRowNumbering()
    {
        return mpRowNumbering;
    }

    inline typename DistributedNumbering<TIndexType>::UniquePointer& pGetColNumbering()
    {
        return mpColNumbering;
    }

    inline typename DistributedVectorImporter<TDataType,TIndexType>::UniquePointer& pGetVectorImporter()
    {
        return mpVectorImporter;
    }

    void SetValue(const TDataType value)
    {
        GetDiagonalBlock().SetValue(value);
        GetOffDiagonalBlock().SetValue(value);
    }

    TIndexType local_size1() const
    {
        return GetDiagonalBlock().size1();
    }

    TIndexType size2() const
    {
        return GetColNumbering().Size();
    }

    inline TIndexType local_nnz() const
    {
        return GetDiagonalBlock().nnz();
    }

    const DataCommunicator& GetComm() const
    {
        return *mpComm;
    }

    const DataCommunicator* pGetComm() const
    {
        return mpComm;
    }

    inline typename CsrMatrix<TDataType,TIndexType>::UniquePointer& pGetDiagonalBlock()
    {
        return mpDiagonalBlock;
    }
    inline typename CsrMatrix<TDataType,TIndexType>::UniquePointer& pGetOffDiagonalBlock()
    {
        return mpOffDiagonalBlock;
    }

    inline CsrMatrix<TDataType,TIndexType>& GetDiagonalBlock()
    {
        return *mpDiagonalBlock;
    }
    inline CsrMatrix<TDataType,TIndexType>& GetOffDiagonalBlock()
    {
        return *mpOffDiagonalBlock;
    }

    inline const CsrMatrix<TDataType,TIndexType>& GetDiagonalBlock() const
    {
        return *mpDiagonalBlock;
    }
    inline const CsrMatrix<TDataType,TIndexType>& GetOffDiagonalBlock() const
    {
        return *mpOffDiagonalBlock;
    }

    const std::map<TIndexType, TIndexType>& GetOffDiagonalLocalIds() const
    {
        return mOffDiagonalLocalIds; //usage: mOffDiagonalLocalIds[global_id] contains the local_id associated to that global_id (for a off diagonal block entry)
    }

    std::map<TIndexType, TIndexType>& GetOffDiagonalLocalIds()
    {
        return mOffDiagonalLocalIds; //usage: mOffDiagonalLocalIds[global_id] contains the local_id associated to that global_id (for a off diagonal block entry)
    }

    const DenseVector<TIndexType>& GetOffDiagonalGlobalIds() const
    {
        return mOffDiagonalGlobalIds;
    }

    DenseVector<TIndexType>& GetOffDiagonalGlobalIds()
    {
        return mOffDiagonalGlobalIds;
    }

    TIndexType GetOffDiagonalBlockLocalId(TIndexType GlobalJ) const
    {
        auto it = mOffDiagonalLocalIds.find(GlobalJ);
        KRATOS_DEBUG_ERROR_IF( it == mOffDiagonalLocalIds.end() ) << "GlobalJ is not in the nonlocal list" << std::endl;
        return it->second;
    }

    TIndexType GetOffDiaGlobalId(TIndexType LocalJ) const
    {
        return mOffDiagonalGlobalIds[LocalJ];
    }

    TDataType& GetLocalDataByGlobalId(TIndexType GlobalI, TIndexType GlobalJ)
    {
        KRATOS_DEBUG_ERROR_IF(  ! GetRowNumbering().IsLocal(GlobalI) ) << "non local row access for GlobalI,GlobalJ = " << GlobalI << " " << GlobalJ << std::endl;

        TIndexType LocalI = GetRowNumbering().LocalId(GlobalI);
        if(GetColNumbering().IsLocal(GlobalJ))
        {
            return GetDiagonalBlock()( LocalI, GetColNumbering().LocalId(GlobalJ) );
        }
        else
        {
            return GetOffDiagonalBlock()( LocalI, GetOffDiagonalBlockLocalId(GlobalJ) );
        }
    }

    TDataType& GetNonLocalCachedDataByGlobalId(TIndexType GlobalI, TIndexType GlobalJ)
    {
        KRATOS_DEBUG_ERROR_IF(  GetRowNumbering().IsLocal(GlobalI) ) << " local row access for GlobalI,GlobalJ = " << GlobalI << " " << GlobalJ << " expected to be nonlocal" << std::endl;
        auto it = mNonLocalData.find(std::make_pair(GlobalI,GlobalJ));
        KRATOS_DEBUG_ERROR_IF(it == mNonLocalData.end()) << " entry GlobalI,GlobalJ = " << GlobalI << " " << GlobalJ << " not found in NonLocalData" << std::endl;
        return it->second;
    }

    DenseVector<TIndexType> GetDiagonalIndex2DataInGlobalNumbering() const
    {
        DenseVector<TIndexType> tmp(GetDiagonalBlock().index2_data().size());
        IndexPartition<TIndexType>(tmp.size()).for_each([&](TIndexType i)
        {
            tmp[i] = GetColNumbering().GlobalId(   GetDiagonalBlock().index2_data()[i]  );
        });
        return tmp;
    }

    DenseVector<TIndexType> GetOffDiagonalIndex2DataInGlobalNumbering() const
    {
        DenseVector<TIndexType> tmp(GetOffDiagonalBlock().index2_data().size());
        IndexPartition<TIndexType>(tmp.size()).for_each([&](TIndexType i)
        {
            tmp[i] = GetOffDiaGlobalId(   GetOffDiagonalBlock().index2_data()[i]  );
        });
        return tmp;
    }

    void ApplyHomogeneousDirichlet(const DistributedSystemVector<TDataType,TIndexType>& rFreeDofsVector,
                                   const TDataType DiagonalValue,
                                   DistributedSystemVector<TDataType,TIndexType>& rRHS)
    {
        KRATOS_TRY
        KRATOS_ERROR_IF(local_size1() != rFreeDofsVector.LocalSize() ) << "ApplyHomogeneousDirichlet: mismatch between row sizes : " << local_size1()
                << " and free_dofs_vector size " << rFreeDofsVector.LocalSize() << std::endl;


        //diagonal block
        GetDiagonalBlock().ApplyHomogeneousDirichlet(rFreeDofsVector.GetLocalData(), DiagonalValue, rRHS.GetLocalData());

        //off diagonal block
        auto off_diag_free = mpVectorImporter->ImportData(rFreeDofsVector); //obtain the nonlocal components of the rFreeDofsVector
        if(GetOffDiagonalBlock().nnz() != 0) {
            KRATOS_ERROR_IF(GetOffDiagonalBlock().size1() != rFreeDofsVector.LocalSize() ) << "ApplyHomogeneousDirichlet - OffDiagonalBlock: mismatch between row sizes : " << GetOffDiagonalBlock().size1()
                    << " and free_dofs_vector size " << rFreeDofsVector.LocalSize() << std::endl;
            KRATOS_ERROR_IF(GetOffDiagonalBlock().size2() != off_diag_free.size() ) << "ApplyHomogeneousDirichlet - OffDiagonalBlock: mismatch between col sizes : " << GetOffDiagonalBlock().size2()
                    << " and free_dofs_vector size " << off_diag_free.size() << std::endl;
            IndexPartition<IndexType>(GetDiagonalBlock().size1()).for_each( [&](IndexType i)
            {
                const IndexType row_begin = GetOffDiagonalBlock().index1_data()[i];
                const IndexType row_end   = GetOffDiagonalBlock().index1_data()[i+1];
                if(std::abs(rFreeDofsVector.GetLocalData()[i] - 1.0) < 1e-14) { //row corresponding to free dofs NOTE that we check if it is approx zero TODO: Epsilon?
                    for(IndexType k = row_begin; k < row_end; ++k) {
                        const IndexType col = GetOffDiagonalBlock().index2_data()[k];
                        GetOffDiagonalBlock().value_data()[k] *= off_diag_free[col]; //note that here we assume that rFreeDofsVector is either 0 or 1
                    }
                } else { //row corresponding to a fixed dof - set to zero all non local row
                    for(IndexType k = row_begin; k < row_end; ++k) {
                        GetOffDiagonalBlock().value_data()[k] = TDataType();
                    }
                }
            });
        }

        KRATOS_CATCH("")
    }
    // TDataType& operator()(TIndexType I, TIndexType J){
    // }

    ///@}
    ///@name Operations
    ///@{
    //TODO
    //+=
    //-=
    //+
    //-
    //*

    // y += A*x  -- where A is *this
    void SpMV(const DistributedSystemVector<TDataType,TIndexType>& x,
              DistributedSystemVector<TDataType,TIndexType>& y) const
    {
        //get off diagonal terms (requires communication)
        auto off_diag_x = mpVectorImporter->ImportData(x);
        GetOffDiagonalBlock().SpMV(off_diag_x,y.GetLocalData());
        GetDiagonalBlock().SpMV(x.GetLocalData(),y.GetLocalData());
    }

    //y = alpha*y + beta*A*x
    void SpMV(const TDataType alpha,
              const DistributedSystemVector<TDataType,TIndexType>& x,
              const TDataType beta,
              DistributedSystemVector<TDataType,TIndexType>& y) const
    {
        //get off diagonal terms (requires communication)
        auto off_diag_x = mpVectorImporter->ImportData(x);
        GetOffDiagonalBlock().SpMV(alpha,off_diag_x,beta,y.GetLocalData());
        GetDiagonalBlock().SpMV(alpha,x.GetLocalData(),beta,y.GetLocalData());
    }

    // y += A^t*x  -- where A is *this
    DistributedVectorExporter<TIndexType>* TransposeSpMV(const DistributedSystemVector<TDataType,TIndexType>& x,
            DistributedSystemVector<TDataType,TIndexType>& y,
            DistributedVectorExporter<TIndexType>* pTransposeExporter = nullptr
                                                        ) const
    {
        GetDiagonalBlock().TransposeSpMV(x.GetLocalData(),y.GetLocalData());

        DenseVector<TDataType> non_local_transpose_data = ZeroVector(GetOffDiagonalBlock().size2());
        GetOffDiagonalBlock().TransposeSpMV(x.GetLocalData(),non_local_transpose_data);

        if(pTransposeExporter == nullptr) //if a nullptr is passed the DistributedVectorExporter is constructed on the flight
        {
            //constructing the exporter requires communication, so the efficiency of the TransposeSpMV can be increased by passing a constructed exporter
            pTransposeExporter = new DistributedVectorExporter<TIndexType>(GetComm(),mOffDiagonalGlobalIds,y.GetNumbering());
        }
        pTransposeExporter->Apply(y,non_local_transpose_data);

        return pTransposeExporter;
    }

    // y += A^t*x  -- where A is *this
    DistributedVectorExporter<TIndexType>* TransposeSpMV(
        TDataType alpha,
        const DistributedSystemVector<TDataType,TIndexType>& x,
        TDataType beta,
        DistributedSystemVector<TDataType,TIndexType>& y,
        DistributedVectorExporter<TIndexType>* pTransposeExporter = nullptr
    ) const
    {
        GetDiagonalBlock().TransposeSpMV(alpha,x.GetLocalData(),beta,y.GetLocalData());

        DenseVector<TDataType> non_local_transpose_data = ZeroVector(GetOffDiagonalBlock().size2());
        GetOffDiagonalBlock().TransposeSpMV(alpha,x.GetLocalData(),beta,non_local_transpose_data);

        if(pTransposeExporter == nullptr) //if a nullptr is passed the DistributedVectorExporter is constructed on the flight
        {
            //constructing the exporter requires communication, so the efficiency of the TransposeSpMV can be increased by passing a constructed exporter
            pTransposeExporter = new DistributedVectorExporter<TIndexType>(GetComm(),mOffDiagonalGlobalIds,y.GetNumbering());
        }
        pTransposeExporter->Apply(y,non_local_transpose_data);

        return pTransposeExporter;
    }

    TDataType NormFrobenius() const
    {
        TDataType diag_norm = GetDiagonalBlock().NormFrobenius();
        TDataType off_diag_norm = GetOffDiagonalBlock().NormFrobenius();
        TDataType sum_squared = std::pow(diag_norm,2) + std::pow(off_diag_norm,2);
        sum_squared = GetComm().SumAll(sum_squared);
        return std::sqrt(sum_squared);
    }


    void BeginAssemble()
    {
        //set to zero non local data
        for(auto& item : mNonLocalData)
            item.second = 0.0;
    }

    void FinalizeAssemble()
    {
        //communicate data to finalize the assembly
        auto& rComm = GetComm();

        std::vector<TDataType> send_data;
        std::vector<TDataType> recv_data;

        //sendrecv data
        for(auto color : mfem_assemble_colors) {
            if(color >= 0) { //-1 would imply no communication
                const auto& direct_senddata_access = mPointersToSendValues[color];
                const auto& direct_recvdata_access = mPointersToRecvValues[color];

                send_data.resize(direct_senddata_access.size());
                recv_data.resize(direct_recvdata_access.size());

                for(TIndexType i=0; i<send_data.size(); ++i) {
                    send_data[i] = *(direct_senddata_access[i]);
                }

                rComm.SendRecv(send_data, color, 0, recv_data, color, 0);

                for(TIndexType i=0; i<recv_data.size(); ++i) {
                    *(direct_recvdata_access[i]) += recv_data[i]; //here we assemble the nonlocal contribution to the local data
                }
            }
        }
    }

    template<class TMatrixType, class TIndexVectorType >
    void Assemble(
        const TMatrixType& rMatrixInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size1() != EquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size2() != EquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;

        for(unsigned int i=0; i<EquationId.size(); ++i) {
            const TIndexType global_i = EquationId[i];
            if(GetRowNumbering().IsLocal(global_i)) {
                for(unsigned int j = 0; j<EquationId.size(); ++j) {
                    const TIndexType global_j = EquationId[j];
                    TDataType& value = GetLocalDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }
            } else {
                for(unsigned int j = 0; j<EquationId.size(); ++j) {
                    const TIndexType global_j = EquationId[j];
                    TDataType& value = GetNonLocalCachedDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }
            }
        }
    }

    void AssembleEntry(const TDataType Value, const TIndexType GlobalI, const TIndexType GlobalJ)
    {
        if(GetRowNumbering().IsLocal(GlobalI)) {
            TDataType& v = GetLocalDataByGlobalId(GlobalI,GlobalJ);
            AtomicAdd(v, Value);
        } else {
            TDataType& v = GetNonLocalCachedDataByGlobalId(GlobalI,GlobalJ);
            AtomicAdd(v, Value);
        }
    }

    template<class TMatrixType, class TIndexVectorType >
    void Assemble(
        const TMatrixType& rMatrixInput,
        const TIndexVectorType& RowEquationId,
        const TIndexVectorType& ColEquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size1() != RowEquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size2() != ColEquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;

        for(unsigned int i=0; i<RowEquationId.size(); ++i) {
            const TIndexType global_i = RowEquationId[i];

            if(GetRowNumbering().IsLocal(global_i)) {
                for(unsigned int j = 0; j<ColEquationId.size(); ++j) {
                    const TIndexType global_j = ColEquationId[j];
                    TDataType& value = GetLocalDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }
            } else {
                for(unsigned int j = 0; j<ColEquationId.size(); ++j) {
                    const TIndexType global_j = ColEquationId[j];
                    TDataType& value = GetNonLocalCachedDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }

            }
        }
    }

    MatrixMapType ToMap() const
    {
        MatrixMapType value_map;
        for(unsigned int i=0; i<local_size1(); ++i) {
            TIndexType row_begin = GetDiagonalBlock().index1_data()[i];
            TIndexType row_end   = GetDiagonalBlock().index1_data()[i+1];
            for(TIndexType k = row_begin; k < row_end; ++k) {
                TIndexType j = GetDiagonalBlock().index2_data()[k];
                TDataType v = GetDiagonalBlock().value_data()[k];
                value_map[ {GetRowNumbering().GlobalId(i),GetColNumbering().GlobalId(j)}] = v;
            }
        }

        for(unsigned int i=0; i<local_size1(); ++i) {
            TIndexType row_begin = GetOffDiagonalBlock().index1_data()[i];
            TIndexType row_end   = GetOffDiagonalBlock().index1_data()[i+1];
            for(TIndexType k = row_begin; k < row_end; ++k) {
                TIndexType j = GetOffDiagonalBlock().index2_data()[k];
                TDataType v = GetOffDiagonalBlock().value_data()[k];
                value_map[ {GetRowNumbering().GlobalId(i),mOffDiagonalGlobalIds[j]}] = v;
            }
        }

        return value_map;
    }

    typename CsrMatrix<TDataType,TIndexType>::Pointer ToSerialCSR(MpiIndexType target_rank=0) const
    {
        // Flatten all data (both indices and values) into a single vector of doubles
        std::vector<double> tmp_data;
        tmp_data.reserve(GetDiagonalBlock().nnz()*3 + GetOffDiagonalBlock().nnz()*3);
        for(unsigned int i=0; i<GetDiagonalBlock().size1(); ++i){
            IndexType row_begin = GetDiagonalBlock().index1_data()[i];
            IndexType row_end   = GetDiagonalBlock().index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                const IndexType j = GetDiagonalBlock().index2_data()[k];
                const TDataType v = GetDiagonalBlock().value_data()[k];
                tmp_data.push_back(GetRowNumbering().GlobalId(i));
                tmp_data.push_back(GetColNumbering().GlobalId(j));
                tmp_data.push_back(v);
            }
        }
        for(unsigned int i=0; i<GetOffDiagonalBlock().size1(); ++i){
            IndexType row_begin = GetOffDiagonalBlock().index1_data()[i];
            IndexType row_end   = GetOffDiagonalBlock().index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                const IndexType j = GetOffDiagonalBlock().index2_data()[k];
                const TDataType v = GetOffDiagonalBlock().value_data()[k];
                tmp_data.push_back(GetRowNumbering().GlobalId(i));
                tmp_data.push_back(mOffDiagonalGlobalIds[j]);
                tmp_data.push_back(v);
            }
        }

        auto collected_data = GetComm().Gatherv(tmp_data,target_rank);

        const MpiIndexType num_processors = GetComm().Size();
        const MpiIndexType my_rank = GetComm().Rank();

        typename SparseContiguousRowGraph<TIndexType>::UniquePointer p_csr_graph;
        typename CsrMatrix<TDataType,TIndexType>::Pointer p_csr_output;

        if(my_rank==target_rank){
            p_csr_graph = std::move( Kratos::make_unique<SparseContiguousRowGraph<TIndexType>>(GetRowNumbering().Size()) );

            for(int i_proc=0; i_proc<num_processors; ++i_proc){
                const auto& data = collected_data[i_proc];
                for(IndexType i=0; i<data.size(); i+=3){
                    IndexType I = static_cast<IndexType>(data[i]);
                    IndexType J = static_cast<IndexType>(data[i+1]);
                    p_csr_graph->AddEntry(I,J);
                }
            }
            p_csr_graph->Finalize();
        }

        if(my_rank==target_rank){
            p_csr_output = Kratos::make_shared<CsrMatrix<TDataType,TIndexType>>(*p_csr_graph);
            p_csr_output->BeginAssemble();
            for(int i_proc=0; i_proc<num_processors; ++i_proc){
                const auto& data = collected_data[i_proc];
                for(IndexType i=0; i<data.size(); i+=3){
                    const IndexType I = static_cast<IndexType>(data[i]);
                    const IndexType J = static_cast<IndexType>(data[i+1]);
                    const double    v = data[i+2];
                    p_csr_output->AssembleEntry(v,I,J);
                }
            }
            p_csr_output->FinalizeAssemble();
        }

        return p_csr_output;
    }

    TDataType NormDiagonal() const
    {
        TDataType diagonal_norm_squared = std::pow(GetDiagonalBlock().NormDiagonal(), 2);
        diagonal_norm_squared = GetComm().SumAll(diagonal_norm_squared);
        return (std::sqrt(diagonal_norm_squared));
    }

    TDataType MaxDiagonal() const
    {
        TDataType diagonal_max = GetDiagonalBlock().MaxDiagonal();
        diagonal_max = GetComm().MaxAll(diagonal_max);
        return diagonal_max;
    }

    TDataType MinDiagonal() const
    {
        TDataType diagonal_min = GetDiagonalBlock().MinDiagonal();
        diagonal_min = GetComm().MinAll(diagonal_min);
        return diagonal_min;
    }

    //TODO
    // LeftScaling
    // RightScaling
    // SymmetricScaling

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DistributedCsrMatrix" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DistributedCsrMatrix" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "--- Diagonal Block: ---" << std::endl;
        rOStream << "size1 : " << GetDiagonalBlock().size1() <<std::endl;
        rOStream << "size2 : " << GetDiagonalBlock().size2() <<std::endl;
        rOStream << "nnz : " << GetDiagonalBlock().nnz() <<std::endl;
        rOStream << "index1_data : " << std::endl;
        for(auto item : GetDiagonalBlock().index1_data())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "index2_data in local numbering: " << std::endl;
        for(auto item : GetDiagonalBlock().index2_data())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "index2_data in global numbering: " << std::endl;
        for(auto item : GetDiagonalIndex2DataInGlobalNumbering())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "value_data  : " << std::endl;
        for(auto item : GetDiagonalBlock().value_data())
            rOStream << item << ",";
        rOStream << std::endl;
        //rOStream << "as map      : " << GetDiagonalBlock().ToMap() << std::endl;
        rOStream << std::endl;
        rOStream << "--- OffDiagonal Block: ---" << std::endl;
        rOStream << "size1 : " << GetOffDiagonalBlock().size1() <<std::endl;
        rOStream << "size2 : " << GetOffDiagonalBlock().size2() <<std::endl;
        rOStream << "nnz : " << GetOffDiagonalBlock().nnz() <<std::endl;
        rOStream << "index1_data : " << std::endl;
        for(auto item : GetOffDiagonalBlock().index1_data())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "index2_data in local numbering: " << std::endl;
        for(auto item : GetOffDiagonalBlock().index2_data())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "index2_data in global numbering: " << std::endl;
        for(auto item : GetOffDiagonalIndex2DataInGlobalNumbering())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "value_data  : " << std::endl;
        for(auto item : GetOffDiagonalBlock().value_data())
            rOStream << item << ",";
        rOStream << std::endl;
        //rOStream << "as map      : " << GetOffDiagonalBlock().ToMap() << std::endl;


    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{
    // A non-recursive binary search function. It returns
    // location of x in given array arr[l..r] is present,
    // otherwise -1
    template< class TVectorType >
    inline TIndexType BinarySearch(const TVectorType& arr,
                                   TIndexType l, TIndexType r, TIndexType x)
    {
        while (l <= r) {
            int m = l + (r - l) / 2;

            // Check if x is present at mid
            if (arr[m] == x)
                return m;

            // If x greater, ignore left half
            if (arr[m] < x)
                l = m + 1;

            // If x is smaller, ignore right half
            else
                r = m - 1;
        }

        // if we reach here, then element was not present
        return -1;
    }

    ///@}
    ///@name Protected Operations
    ///@{
    void PrepareNonLocalCommunications(const DistributedSparseGraph<TIndexType>& rSparseGraph)
    {
        auto& rComm = GetComm();
        const auto& nonlocal_graphs = rSparseGraph.GetNonLocalGraphs();
        std::vector<int> send_list;
        for(unsigned int id = 0; id<nonlocal_graphs.size(); ++id)
            if( !nonlocal_graphs[id].IsEmpty())
                send_list.push_back(id);

        mfem_assemble_colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, rComm);

        //sendrecv data
        for(auto color : mfem_assemble_colors) {
            if(color >= 0) { //-1 would imply no communication
                const auto& send_graph = nonlocal_graphs[color];
                auto& direct_senddata_access = mPointersToSendValues[color];
                auto& send_ij = mSendCachedIJ[color];

                for(auto row_it=send_graph.begin(); row_it!=send_graph.end(); ++row_it) {
                    const auto remote_local_I = row_it.GetRowIndex();
                    const auto remote_global_I = GetRowNumbering().RemoteGlobalId(remote_local_I, color);

                    for(auto J : *row_it) {
                        TDataType& value = mNonLocalData[std::make_pair(remote_global_I,J)]; //here we create the I,J entry in the nonlocal data (entry was there in the graph!)
                        direct_senddata_access.push_back(&value); //storing a direct pointer to the value contained in the data structure
                        send_ij.push_back(remote_global_I);
                        send_ij.push_back(J);
                    }
                }

                //NOTE: this can be made nonblocking
                mRecvCachedIJ[color] = rComm.SendRecv(send_ij, color, color);
                auto& recv_ij = mRecvCachedIJ[color];

                auto& direct_recvdata_access = mPointersToRecvValues[color];

                for(TIndexType k=0; k<recv_ij.size(); k+=2) {
                    TIndexType I = recv_ij[k];
                    TIndexType J = recv_ij[k+1];
                    auto& value = GetLocalDataByGlobalId(I,J);
                    direct_recvdata_access.push_back(&value);
                }
            }
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    const DataCommunicator* mpComm;

    typename DistributedNumbering<TIndexType>::UniquePointer mpRowNumbering;
    typename DistributedNumbering<TIndexType>::UniquePointer mpColNumbering;

    typename BlockMatrixType::UniquePointer mpDiagonalBlock = Kratos::make_unique< BlockMatrixType >();
    typename BlockMatrixType::UniquePointer mpOffDiagonalBlock = Kratos::make_unique< BlockMatrixType >();
    MatrixMapType mNonLocalData; //data which is assembled locally and needs to be communicated to the owner

    //this map tells for an index J which does not belong to the local diagonal block which is the corresponding localJ
    std::map<TIndexType, TIndexType> mOffDiagonalLocalIds; //usage: mOffDiagonalLocalIds[global_id] contains the local_id associated to that global_id (for a off diagonal block entry)
    DenseVector<TIndexType> mOffDiagonalGlobalIds; //usage: mOffDiagonalGlobalIds[local_id] contains the global_id associated

    std::vector<int> mfem_assemble_colors; //coloring of communication
    std::unordered_map< unsigned int, std::vector<TIndexType> > mRecvCachedIJ; //recv_ij contains i,j to receive one after the other
    std::unordered_map< unsigned int, std::vector<TIndexType> > mSendCachedIJ; //recv_ij contains i,j to receive one after the other
    std::unordered_map< unsigned int, std::vector<TDataType*> > mPointersToRecvValues; //this contains direct pointers into the data contained in mNonLocalData, prepared so to speed up communications
    std::unordered_map< unsigned int, std::vector<TDataType*> > mPointersToSendValues; //this contains direct pointers into mDiagonalBlock and mOffDiagonalBlock, prepared so to speed up communication

    std::unique_ptr<DistributedVectorImporter<TDataType,TIndexType>> mpVectorImporter;

    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("CommunicatorName",ParallelEnvironment::RetrieveRegisteredName(*mpComm));
        rSerializer.save("RowNumbering",mpRowNumbering);
        rSerializer.save("ColNumbering",mpColNumbering);
        rSerializer.save("mpDiagonalBlock",mpDiagonalBlock);
        rSerializer.save("mpOffDiagonalBlock",mpOffDiagonalBlock);
        rSerializer.save("mNonLocalData",mNonLocalData);
        rSerializer.save("mSendCachedIJ",mSendCachedIJ);
        rSerializer.save("mRecvCachedIJ",mRecvCachedIJ);
        rSerializer.save("mOffDiagonalLocalIds",mOffDiagonalLocalIds);
        rSerializer.save("mfem_assemble_colors",mfem_assemble_colors);
        rSerializer.save("mpVectorImporter",mpVectorImporter);
        //note that the direct access vectors are not saved, we will need to reconstruct them basing on send_ij and recv_ij

    }

    void load(Serializer& rSerializer)
    {
        std::string comm_name;
        rSerializer.load("CommunicatorName",comm_name);
        mpComm = &ParallelEnvironment::GetDataCommunicator(comm_name);
        rSerializer.load("RowNumbering",mpRowNumbering);
        rSerializer.load("ColNumbering",mpColNumbering);
        rSerializer.load("mpDiagonalBlock",mpDiagonalBlock);
        rSerializer.load("mpOffDiagonalBlock",mpOffDiagonalBlock);
        rSerializer.load("mNonLocalData",mNonLocalData);
        rSerializer.load("mSendCachedIJ",mSendCachedIJ);
        rSerializer.load("mRecvCachedIJ",mRecvCachedIJ);
        rSerializer.load("mOffDiagonalLocalIds",mOffDiagonalLocalIds);
        rSerializer.load("mfem_assemble_colors",mfem_assemble_colors);
        rSerializer.load("mpVectorImporter",mpVectorImporter);

        ReconstructDirectAccessVectors();
    }


    ///@}
    ///@name Private Operations
    ///@{
    void ReconstructDirectAccessVectors()
    {
        //compute direct pointers to data
        for(auto color : mfem_assemble_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                const auto& send_ij = mSendCachedIJ[color];
                const auto& recv_ij = mRecvCachedIJ[color];

                auto& direct_senddata_access = mPointersToRecvValues[color];
                for(TIndexType i=0; i<send_ij.size(); i+=2)
                {
                    TIndexType I = recv_ij[i];
                    TIndexType J = recv_ij[i+1];
                    auto& value = GetNonLocalCachedDataByGlobalId(I,J);
                    direct_senddata_access.push_back(&value);
                }

                auto& direct_recvdata_access = mPointersToRecvValues[color];
                for(TIndexType k=0; k<recv_ij.size(); k+=2)
                {
                    TIndexType I = recv_ij[k];
                    TIndexType J = recv_ij[k+1];
                    auto& value = GetLocalDataByGlobalId(I,J);
                    direct_recvdata_access.push_back(&value);
                }
            }
        }
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{



    ///@}

}; // Class DistributedCsrMatrix

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDataType>
inline std::istream& operator >> (std::istream& rIStream,
                                  DistributedCsrMatrix<TDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DistributedCsrMatrix<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_CSR_MATRIX_H_INCLUDED  defined


