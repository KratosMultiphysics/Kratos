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
#if !defined(KRATOS_SYSTEM_VECTOR_H_INCLUDED )
#define  KRATOS_SYSTEM_VECTOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_utilities.h"


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

/// Provides a SystemVector which implements FEM assemble capabilities, as well as some vector operations
template<class TDataType=double, class TIndexType=std::size_t>
class SystemVector final
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;

    /// Pointer definition of SystemVector
    KRATOS_CLASS_POINTER_DEFINITION(SystemVector);

    ///@}
    ///@name Life Cycle
    ///@{

    SystemVector(const SparseGraph<IndexType>& rGraph){
        mpComm = rGraph.pGetComm();
        mData.resize(rGraph.Size(),false);
    }

    SystemVector(const SparseContiguousRowGraph<IndexType>& rGraph){
        mpComm = rGraph.pGetComm();
        mData.resize(rGraph.Size(),false);
    }

    SystemVector(IndexType size, DataCommunicator& rComm=ParallelEnvironment::GetDataCommunicator("Serial")){
        if(rComm.IsDistributed())
            KRATOS_ERROR << "Attempting to construct a serial system_vector with a distributed communicator" << std::endl;
        mpComm = &rComm;
        mData.resize(size,false);
    }

    SystemVector(
        const Vector& data,
        DataCommunicator& rComm = ParallelEnvironment::GetDataCommunicator("Serial")) {
        if(rComm.IsDistributed())
            KRATOS_ERROR << "Attempting to construct a serial system_vector with a distributed communicator" << std::endl;
        mpComm = &rComm;
        mData.resize(data.size(),false);
        noalias(mData) = data;
    }
    /// Copy constructor.
    explicit SystemVector(const SystemVector<TDataType,TIndexType>& rOtherVector){
        mpComm = rOtherVector.mpComm;
        mData.resize(rOtherVector.size(),false);

        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] = rOtherVector[i];
        });
    }

    /// Destructor.
    ~SystemVector(){}

    const DataCommunicator& GetComm() const
    {
        return *mpComm;
    }

    const DataCommunicator* pGetComm() const
    {
        return mpComm;
    }

    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
        mData.clear();
    }

    void SetValue(const TDataType value)
    {
        IndexPartition<IndexType>(mData.size()).for_each([&](IndexType i){
            mData[i] = value;
        });
    }

    IndexType size() const
    {
        return mData.size();
    }

    TDataType& operator()(IndexType I){
        return mData[I];
    }

    const TDataType& operator()(IndexType I) const{
        return mData[I];
    }

    TDataType& operator[](IndexType I){
        return mData[I];
    }

    const TDataType& operator[](IndexType I) const{
        return mData[I];
    }

    ///provides low level access to internal data
    DenseVector<TDataType>& data()
    {
        return mData;
    }

    ///provides low level access to internal data
    const DenseVector<TDataType>& data() const
    {
        return mData;
    }

    void Add(const TDataType factor,
             const SystemVector& rOtherVector
            )
    {
        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] += factor*rOtherVector[i];
        });
    }

    /// Assignment operator.
    SystemVector& operator=(SystemVector const& rOtherVector){
        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] = rOtherVector[i];
        });
        return *this;
    }


    SystemVector& operator+=(const SystemVector& rOtherVector)
    {
        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] += rOtherVector[i];
        });
        return *this;
    }

    SystemVector& operator-=(const SystemVector& rOtherVector)
    {
        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] -= rOtherVector[i];
        });
        return *this;
    }

    SystemVector& operator*=(const TDataType multiplier_factor)
    {
        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] *= multiplier_factor;
        });
        return *this;
    }

    SystemVector& operator/=(const TDataType divide_factor)
    {
        IndexPartition<IndexType>(size()).for_each([&](IndexType i){
            (*this)[i] /= divide_factor;
        });
        return *this;
    }

    TDataType Dot(const SystemVector& rOtherVector, IndexType gather_on_rank=0)
    {
        KRATOS_WARNING_IF("SystemVector", gather_on_rank != 0) << "the parameter gather_on_rank essentially does nothing for a non-distribued vector. It is added to have the same interface as for the distributed_system_vector" << std::endl;

        auto partition = IndexPartition<IndexType>(size());
        TDataType dot_value = partition.template for_each< SumReduction<TDataType> >([&](IndexType i){
                return (*this)[i]*rOtherVector[i];
            });

        return dot_value; // note that the value to be reduced should be returned
    }


    ///@}
    ///@name Operations
    ///@{
    void BeginAssemble(){} //the SMP version does nothing. This function is there to be implemented in the MPI case

    void FinalizeAssemble(){} //the SMP version does nothing. This function is there to be implemented in the MPI case

    template<class TVectorType, class TIndexVectorType >
    void Assemble(
        const TVectorType& rVectorInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rVectorInput.size() != EquationId.size());

        for(unsigned int i=0; i<EquationId.size(); ++i){
            IndexType global_i = EquationId[i];
            AssembleEntry(rVectorInput[i], global_i);
        }
    }

    void AssembleEntry(
        const TDataType rValue,
        const IndexType GlobalI)
    {
        KRATOS_DEBUG_ERROR_IF(GlobalI > mData.size());
        AtomicAdd(mData[GlobalI] , rValue);
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
        buffer << "SystemVector" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << "SystemVector" << std::endl;
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        std::cout << mData << std::endl;
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
    DenseVector<TDataType> mData;

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

}; // Class SystemVector

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType, class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                SystemVector<TDataType,TIndexType>& rThis)
                {
                    return rIStream;
                }

/// output stream function
template<class TDataType, class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                const SystemVector<TDataType,TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SYSTEM_VECTOR_H_INCLUDED  defined


