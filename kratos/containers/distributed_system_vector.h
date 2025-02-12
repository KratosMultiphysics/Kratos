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
#if !defined(KRATOS_DISTRIBUTED_SYSTEM_VECTOR_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_SYSTEM_VECTOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <functional>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/distributed_system_vector.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_numbering.h"
#include "containers/distributed_vector_exporter.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Provides a DistributedSystemVector which implements FEM assemble capabilities
template<class TDataType=double, class TIndexType=std::size_t>
class DistributedSystemVector final
{
public:
    ///@name Type Definitions
    ///@{
    typedef TDataType value_type;
    typedef TIndexType IndexType;
    typedef int MpiIndexType;


    /// Pointer definition of DistributedSystemVector
    KRATOS_CLASS_POINTER_DEFINITION(DistributedSystemVector);

    ///@}
    ///@name Life Cycle
    ///@{

    DistributedSystemVector(const DistributedSparseGraph<IndexType>& rGraph)
            :
            mpComm(&rGraph.GetComm())
    {
        mpNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >( rGraph.GetRowNumbering());

        mLocalData.resize(rGraph.LocalSize(),false);

        //now add entries in nonlocal data;
        const auto& r_non_local_graphs = rGraph.GetNonLocalGraphs();

        for(IndexType cpu_id=0; cpu_id < r_non_local_graphs.size(); cpu_id++) //this loop cannot be done in parallel since we need to allocate memory in the unordered_map
        {
            const auto& graph = r_non_local_graphs[cpu_id];
            for(auto item=graph.begin(); item!=graph.end(); ++item)
            {
                    IndexType global_id = GetNumbering().RemoteGlobalId(item.GetRowIndex(), cpu_id);
                    mNonLocalData[global_id] = TDataType(); //first touching of nonlocaldata
            }
        }
    }

    /// Copy constructor.
    explicit DistributedSystemVector(DistributedSystemVector const& rOther)
    :
    mpComm(&rOther.GetComm())
    {
        mpNumbering = Kratos::make_unique<DistributedNumbering<IndexType>>(rOther.GetNumbering());
        KRATOS_ERROR_IF(LocalSize() != rOther.LocalSize());
        KRATOS_ERROR_IF(Size() != rOther.Size());
        //copying the data
        mLocalData.resize(rOther.LocalSize(),false);
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] = rOther[i];
        });
        mNonLocalData = rOther.mNonLocalData;
        if(rOther.mpexporter != nullptr) //this will happen if Finalize was not called on the vector we construct from
            mpexporter = Kratos::make_unique<DistributedVectorExporter<IndexType>>(rOther.GetExporter());
    }

    //WARNING: if a Distributed vector is constructed using this constructor, it cannot be used for assembly (non local entries are not precomputed)
    DistributedSystemVector(const DistributedNumbering<IndexType>& rNumbering)
            :
            mpComm(&rNumbering.GetComm())
    {
        mpNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >( rNumbering );
        mLocalData.resize(rNumbering.LocalSize(),false);
    }

    /// Destructor.
    ~DistributedSystemVector(){}

    ///@}
    ///@name Operators
    ///@{
    const DataCommunicator& GetComm() const{
        return *mpComm;
    }

    const DataCommunicator* pGetComm() const{
        return mpComm;
    }

    const DistributedNumbering<IndexType>& GetNumbering() const
    {
        return *mpNumbering;
    }

    void Clear()
    {
        mLocalData.clear();
    }

    void SetValue(const TDataType value)
    {
        IndexPartition<IndexType>(mLocalData.size()).for_each([&](IndexType i){
            mLocalData[i] = value;
        });
    }

    TDataType& operator()(IndexType I){
        return mLocalData[I];
    }


    const TDataType& operator()(IndexType I) const{
        return mLocalData[I];
    }

    TDataType& operator[](IndexType I){
        return mLocalData[I];
    }

    const TDataType& operator[](IndexType I) const{
        return mLocalData[I];
    }

    inline IndexType Size() const{
        return mpNumbering->Size();
    }

    inline IndexType LocalSize() const{
        return mpNumbering->LocalSize();
    }

    DenseVector<TDataType>& GetLocalData(){
        return mLocalData;
    }

    const DenseVector<TDataType>& GetLocalData() const{
        return mLocalData;
    }

    const DistributedVectorExporter<TIndexType>& GetExporter() const{
        KRATOS_DEBUG_ERROR_IF(mpexporter==nullptr) << " mpexporter was not initialized, GetExporter() cannot be used" << std::endl;
        return *mpexporter;
    }

    //function to add and empty entry to mNonLocalData
    void AddEntry(TIndexType GlobalI) //WARNING: NOT THREADSAFE!
    {
        if( !GetNumbering().IsLocal(GlobalI))
            mNonLocalData[GlobalI] = TDataType();
    }

    //function to add empty entries entry to mNonLocalData
    template<class TIteratorType>
    void AddEntries(TIteratorType it_begin, TIteratorType it_end) //WARNING: NOT THREADSAFE!
    {
        for(TIteratorType it=it_begin; it!=it_end; ++it)
            AddEntry(*it);
    }

    TDataType Norm() const
    {
        TDataType norm_squared = IndexPartition<TIndexType>(LocalSize()).template for_each< SumReduction<TDataType> >( [this](TIndexType i){
            return std::pow((mLocalData)[i], 2);
        });
        norm_squared = GetComm().SumAll(norm_squared);
        return std::sqrt(norm_squared);
    }

    TDataType Dot(const DistributedSystemVector& rOtherVector, MpiIndexType gather_on_rank=0)
    {
        const auto& other_data = rOtherVector.GetLocalData();
        TDataType dot_value = IndexPartition<IndexType>(mLocalData.size()).template for_each<SumReduction<TDataType>>([&](IndexType i){
                return mLocalData[i]*other_data[i];
            });

        dot_value = GetComm().Sum(dot_value, gather_on_rank);
        if(GetComm().Rank() != gather_on_rank) dot_value = -1; //give an impossible result in case it is not on the reduction rank
        return dot_value; // note that the value to be reduced should be returned
    }


    void Add(const TDataType factor,
             const DistributedSystemVector& rOtherVector
            )
    {
        KRATOS_ERROR_IF(LocalSize() != rOtherVector.LocalSize()) << "size mismatch in Add Function. LocalSize is " << LocalSize()
                << " " << " rOtherVector.LocalSize() " << rOtherVector.LocalSize() << std::endl;

        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] += factor*rOtherVector.mLocalData[i];
        });
    }

    /// Assignment operator.
    DistributedSystemVector& operator=(DistributedSystemVector const& rOtherVector){
        KRATOS_ERROR_IF(LocalSize() != rOtherVector.LocalSize()) << "size mismatch in assinement operator. LocalSize is " << LocalSize()
                << " " << " rOtherVector.LocalSize() " << rOtherVector.LocalSize() << std::endl;
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] = rOtherVector.mLocalData[i];
        });
        return *this;
    }


    DistributedSystemVector& operator+=(const DistributedSystemVector& rOtherVector)
    {
        KRATOS_ERROR_IF(LocalSize() != rOtherVector.LocalSize()) << "size mismatch in += operator. LocalSize is " << LocalSize()
                << " " << " rOtherVector.LocalSize() " << rOtherVector.LocalSize() << std::endl;
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] += rOtherVector.mLocalData[i];
        });
        return *this;
    }

    DistributedSystemVector& operator-=(const DistributedSystemVector& rOtherVector)
    {
        KRATOS_ERROR_IF(LocalSize() != rOtherVector.LocalSize()) << "size mismatch in -= operator. LocalSize is " << LocalSize()
                << " " << " rOtherVector.LocalSize() " << rOtherVector.LocalSize() << std::endl;
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] -= rOtherVector.mLocalData[i];
        });
        return *this;
    }

    DistributedSystemVector& operator*=(const TDataType& multiplier_factor)
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] *= multiplier_factor;
        });
        return *this;
    }

    DistributedSystemVector& operator/=(const TDataType& divide_factor)
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] /= divide_factor;
        });
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
    void BeginAssemble(){

        //mount exporter. After this point it will be all threadsafe since we will consider as frozen the mNonLocalData structure
        DenseVector<double> non_local_global_ids(mNonLocalData.size());
        IndexType counter = 0;
        for(auto & item : mNonLocalData) //cannot be done in parallel
        {
            non_local_global_ids[counter] = item.first;
            item.second = 0.0; //set to zero prior to assmeply
            counter++;
        }
        auto ptmp = Kratos::make_unique< DistributedVectorExporter<TIndexType> >(GetComm(), non_local_global_ids, GetNumbering());
        mpexporter.swap(ptmp);
    }

    //communicate data to finalize the assembly
    void FinalizeAssemble()
    {
        DenseVector<double> non_local_data(mNonLocalData.size());

        IndexType counter = 0;
        for(auto & item : mNonLocalData) //cannot be done in parallel
        {
            non_local_data[counter] = item.second;
            counter++;
        }

        mpexporter->Apply(*this,non_local_data);
    }

    template<class TVectorType, class TIndexVectorType >
    void Assemble(
        const TVectorType& rVectorInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rVectorInput.size() != EquationId.size());

        for(unsigned int i=0; i<EquationId.size(); ++i){
            IndexType global_i = EquationId[i];

            if(GetNumbering().IsLocal(global_i))
            {
                IndexType local_i = GetNumbering().LocalId(global_i);
                AtomicAdd(mLocalData(local_i) , rVectorInput[i]);
            }
            else
            {
                auto it = (mNonLocalData.find( global_i ));
                KRATOS_DEBUG_ERROR_IF(it == mNonLocalData.end()) << "global_i = "<< global_i << " not in mNonLocalData. Please note that the mNonLocalData structure is assumed fixed after InitializeAssemble is called" << std::endl; //this error is needed for threadsafety
                TDataType& value = (*it).second;
                AtomicAdd(value , rVectorInput[i]);
            }
        }
    }

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
        buffer << "DistributedSystemVector" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DistributedSystemVector LOCAL DATA:" << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        const auto k = GetComm().Rank();
        std::cout << "Local Size=" << LocalSize() << " Total Size=" << GetNumbering().Size() << std::endl;
        std::cout << "Numbering begins with " << GetNumbering().GetCpuBounds()[k] << " ends at " << GetNumbering().GetCpuBounds()[k+1] << std::endl;
        std::cout << mLocalData << std::endl;
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


    ///@}
    ///@name Protected Operations
    ///@{


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
    typename DistributedNumbering<IndexType>::UniquePointer mpNumbering;

    DenseVector<TDataType> mLocalData; //contains the local data
    std::unordered_map<IndexType, TDataType> mNonLocalData; //contains non local data as {global_id,value}

    typename DistributedVectorExporter<TIndexType>::UniquePointer mpexporter = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class DistributedSystemVector

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType, class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                DistributedSystemVector<TDataType,TIndexType>& rThis)
                {
                    return rIStream;
                }

/// output stream function
template<class TDataType, class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                const DistributedSystemVector<TDataType,TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_SYSTEM_VECTOR_H_INCLUDED  defined


