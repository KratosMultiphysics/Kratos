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
#if !defined(KRATOS_DISTRIBUTED_NUMBERING_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_NUMBERING_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
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

///This function provides essential capabilities for mapping between local and global ids in a distributed vector
template<class TIndexType=std::size_t>
class DistributedNumbering
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef int MpiIndexType;

    /// Pointer definition of DistributedNumbering
    KRATOS_CLASS_POINTER_DEFINITION(DistributedNumbering);

    ///@}
    ///@name Life Cycle
    ///@{

    //constructor by local size - computes cpu bounds automatically
    DistributedNumbering(
        const DataCommunicator& rComm,
        const IndexType LocalSize
        )
        : DistributedNumbering(&rComm,LocalSize)
    {}

    DistributedNumbering(
        const DataCommunicator* pComm,
        const IndexType LocalSize)
        :
        mpComm(pComm)
    {
        mCpuBounds.resize(pComm->Size()+1);

        std::vector<IndexType> send_vect{LocalSize};
        std::vector<IndexType> local_sizes = pComm->AllGather(send_vect);
        mCpuBounds[0] = 0;
        for(unsigned int i=1; i<mCpuBounds.size(); ++i)
            mCpuBounds[i] = mCpuBounds[i-1] + local_sizes[i-1];
    }

    //constructor by global sizes and nranks - Requires no communication
    DistributedNumbering(
        const DataCommunicator& rComm,
        const IndexType TotalSize,
        const MpiIndexType Nranks //note that we could obtain this by the rComm, but we need to distinguish this constructor from the previous
        )
        : DistributedNumbering(&rComm, TotalSize,Nranks)
    {}

    DistributedNumbering(
        const DataCommunicator* pComm,
        const IndexType TotalSize,
        const MpiIndexType Nranks //note that we could obtain this by the rComm, but we need to distinguish this constructor from the previous
        )
        :
        mpComm(pComm)
    {
        KRATOS_ERROR_IF(pComm->Size() != Nranks) << "We expect Nranks to be the same as pComm->size()" << std::endl;

        mCpuBounds.resize(pComm->Size()+1);

        const IndexType local_size = TotalSize / Nranks;
        mCpuBounds[0] = 0;
        mCpuBounds[Nranks] = TotalSize;
        for (int i=1; i<Nranks; i++) {
            mCpuBounds[i] = mCpuBounds[i-1] + local_size;
        }
    }

    DistributedNumbering(
        const DataCommunicator* pComm,
        const std::vector<IndexType>& CpuBounds)
        :
        mpComm(pComm), mCpuBounds(CpuBounds)
    {
    }

    DistributedNumbering(
        const DataCommunicator& rComm,
        const std::vector<IndexType>& CpuBounds)
        :
        mpComm(&rComm), mCpuBounds(CpuBounds)
    {
    }

    ///Copy constructor
    DistributedNumbering(const DistributedNumbering& rOther)
        :
        mpComm(rOther.mpComm), mCpuBounds(rOther.GetCpuBounds())
    {
    }

    /// Destructor.
    virtual ~DistributedNumbering() {}

    ///@}
    ///@name Operators
    ///@{
    const DataCommunicator& GetComm() const
    {
        return *mpComm;
    }

    const DataCommunicator* pGetComm() const
    {
        return mpComm;
    }

    inline IndexType LocalSize() const
    {
        const auto k = GetComm().Rank();
        return mCpuBounds[k+1]-mCpuBounds[k];
    }

    inline IndexType Size() const
    {
        return mCpuBounds[mCpuBounds.size()-1];
    }

    inline bool IsLocal(const IndexType I) const
    {
        const auto k = GetComm().Rank();
        return (I>=mCpuBounds[k] && I<mCpuBounds[k+1]);
    }

    inline IndexType LocalId(const IndexType rGlobalId) const
    {
        const auto k = GetComm().Rank();
        return rGlobalId-mCpuBounds[k];
    }

    inline IndexType GlobalId(const IndexType rLocalId) const
    {
        const auto k = GetComm().Rank();
        return rLocalId+mCpuBounds[k];
    }

    inline IndexType RemoteLocalId(const IndexType rGlobalId, const IndexType rOwnerRank) const
    {
        return rGlobalId-mCpuBounds[rOwnerRank];
    }

    inline IndexType RemoteGlobalId(const IndexType rRemoteLocalId, const IndexType rOwnerRank) const
    {
        return rRemoteLocalId+mCpuBounds[rOwnerRank];
    }

    inline IndexType OwnerRank(const IndexType RowIndex) const
    {
        //position of element just larger than RowIndex in mCpuBounds
        auto it = std::upper_bound(mCpuBounds.begin(), mCpuBounds.end(), RowIndex);

        KRATOS_DEBUG_ERROR_IF(it == mCpuBounds.end()) <<
                "row RowIndex " << RowIndex <<
                " is not owned by any processor " << std::endl;

        IndexType owner_rank = (it-mCpuBounds.begin()-1);

        return owner_rank;

    }

    const std::vector<IndexType>& GetCpuBounds() const
    {
        return mCpuBounds;
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DistributedNumbering: CpuBounds = " << mCpuBounds << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DistributedNumbering";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


protected:

private:
    const DataCommunicator* mpComm;
    std::vector<IndexType> mCpuBounds;

    /// Assignment operator.
    DistributedNumbering& operator=(DistributedNumbering const& rOther)
    {
        mpComm = rOther.mpComm;
        this->GetCpuBounds() = rOther.GetCpuBounds();
    }

    ///@}

}; // Class DistributedNumbering

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType, class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                                  DistributedNumbering<TIndexType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TDataType, class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DistributedNumbering<TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_NUMBERING_H_INCLUDED  defined


