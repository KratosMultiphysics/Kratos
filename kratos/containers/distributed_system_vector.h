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


// External includes


// Project includes
#include "includes/define.h"
#include "containers/distributed_system_vector.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

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
class DistributedSystemVector
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;

    /// Pointer definition of DistributedSystemVector
    KRATOS_CLASS_POINTER_DEFINITION(DistributedSystemVector);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistributedSystemVector(const DistributedSparseGraph& rGraph)
            :
            mrComm(rGraph.GetComm()),
            mLocalBounds(rGraph.GetLocalBounds()),
            mCpuBounds(rGraph.GetCpuBounds())
    {
        mLocalData.resize(rGraph.Size(),false);


        //now add entries in nonlocal data;
        const auto& r_non_local_graphs = rGraph.GetNonLocalGraphs();
        for(unsigned int cpu_id = 0; cpu_id<r_non_local_graphs.size(); ++cpu_id)
        IndexPartition<IndexType>(r_non_local_graphs.size()).for_each([&](IndexType cpu_id)
        {
            for(const auto& item : r_non_local_graphs[cpu_id])
            {
                IndexType row = item.first;
                mNonLocalData[cpu_id][row] = TDataType(); //first touching of nonlocaldata
            }
        });
    }

    /// Destructor.
    virtual ~DistributedSystemVector(){}

    ///@}
    ///@name Operators
    ///@{
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

    IndexType Size() const
    {
        return mLocalData.size();
    }

    TDataType& operator()(IndexType I){
        return mLocalData[I];
    }


    const TDataType& operator()(IndexType I) const{
        return mLocalData[I];
    }

    bool IsLocal(const IndexType I) const
    {
        return (I>=mLocalBounds[0] && I<mLocalBounds[1]);
    }

    IndexType LocalId(const IndexType rGlobalId) const
    {
        return rGlobalId-mLocalBounds[0];
    }

    IndexType GlobalId(const IndexType rLocalId) const
    {
        return rLocalId+mLocalBounds[0];
    }

    IndexType RemoteLocalId(const IndexType rGlobalId, const IndexType rOwnerRank) const
    {
        return rGlobalId-mCpuBounds[rOwnerRank];
    } 

    IndexType RemoteGlobalId(const IndexType rRemoteLocalId, const IndexType rOwnerRank) const
    {
        return rRemoteLocalId+mCpuBounds[rOwnerRank];
    } 

    IndexType OwnerRank(const IndexType RowIndex)
    {
        //position of element just larger than RowIndex in mCpuBounds
        auto it = std::upper_bound(mCpuBounds.begin(), mCpuBounds.end(), RowIndex);

        KRATOS_DEBUG_ERROR_IF(it == mCpuBounds.end()) <<
            "row RowIndex " << RowIndex <<
            " is not owned by any processor " << std::endl;

        IndexType owner_rank = (it-mCpuBounds.begin()-1);

        KRATOS_DEBUG_ERROR_IF(owner_rank < 0) <<
            "row RowIndex " << RowIndex <<
            " is not owned by any processor " << std::endl;

        return owner_rank;

    }


    ///@}
    ///@name Operations
    ///@{
    void BeginAssemble(){
        //TODO set to zero nonlocal data prior to assembly
    } 

    void FinalizeAssemble(){
        //TODO communications
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
            
            if(Islocal(global_i))
            {
                IndexType local_i = LocalId(global_i);
                AtomicAdd(mLocalData(local_i) , rVectorInput[i]);
            }
            else
            {
                auto owner_rank = OwnerRank(global_i);
                auto& it = *(mNonLocalData[owner_rank].find(global_i));
                KRATOS_DEBUG_ERROR_IF(it == mNonLocalData[owner_rank].end()) << "global_i = "<< global_i << " not in mNonLocalData" << std::endl;
                TDataType& value = *it;
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
    virtual std::string Info() const
    {
std::stringstream buffer;
    buffer << "DistributedSystemVector" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DistributedSystemVector";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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
    const DataCommunicator& mrComm;
    std::vector<IndexType> mLocalBounds;
    std::vector<IndexType> mCpuBounds;

    DenseVector<TDataType> mLocalData; //contains the local data
    std::unordered_map<IndexType, TDataType> mNonLocalData;

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

    /// Assignment operator.
    DistributedSystemVector& operator=(DistributedSystemVector const& rOther){}

    /// Copy constructor.
    DistributedSystemVector(DistributedSystemVector const& rOther){}

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


