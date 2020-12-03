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
#if !defined(KRATOS_DISTRIBUTED_VECTOR_EXPORTER_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_VECTOR_EXPORTER_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <functional>


// External includes


// Project includes
#include "includes/define.h"
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

/// Provides a DistributedVectorExporter which implements FEM assemble capabilities
template<class TIndexType=std::size_t>
class DistributedVectorExporter
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef int MpiIndexType;

    /// Pointer definition of DistributedVectorExporter
    KRATOS_CLASS_POINTER_DEFINITION(DistributedVectorExporter);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    template< class TGlobalIndicesVectorType>
    DistributedVectorExporter(
        const DataCommunicator& rComm,
        const TGlobalIndicesVectorType& rGlobalIndices,
        const DistributedNumbering<IndexType>& rNumbering)
        :
        mrComm(rComm)
    {
        mpNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >(rNumbering);
        std::unordered_map<int, std::vector<IndexType>> to_send_remote_local_id; //do not need to store this
        for(unsigned int local_i=0; local_i<rGlobalIndices.size(); ++local_i)
        {
            IndexType global_i = rGlobalIndices[local_i];
            MpiIndexType owner_rank = mpNumbering->OwnerRank(global_i);
            IndexType remote_local_i = mpNumbering->RemoteLocalId(global_i, owner_rank);
            mPositionWithinData[owner_rank].push_back(local_i);
            to_send_remote_local_id[owner_rank].push_back(remote_local_i);
        }
        mToRecvLocalId[GetComm().Rank()] = std::move(to_send_remote_local_id[GetComm().Rank()]); //for the current rank do a memcopy

        

        //compute communication plan
        std::vector<MpiIndexType> send_list;
        for(const auto& item : to_send_remote_local_id)
        {
            MpiIndexType cpu_id = item.first;
            if(cpu_id != GetComm().Rank())
                send_list.push_back(cpu_id);
        }
        mVectorCommColors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, rComm);

        //communicate the remote_local_id so that the other node knows what to send
        for(auto color : mVectorCommColors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                //NOTE: this can be made nonblocking 
                mToRecvLocalId[color] = rComm.SendRecv(to_send_remote_local_id[color], color, color); //TODO, we know all the sizes, we shall use that!
            }
        }

        //ensure that entries exist for each color involved
        for(auto color : mVectorCommColors)
        {
            if(color >= 0)
            {
                //create if they do not exist
                mToRecvLocalId[color]; 
                mPositionWithinData[color];
            }
        }
        mToRecvLocalId[GetComm().Rank()]; //create if it does not exist
        mPositionWithinData[GetComm().Rank()]; //create if it does not exist
    }

    /// Copy constructor.
    DistributedVectorExporter(DistributedVectorExporter const& rOther)
        : 
        mrComm(rOther.mrComm),
        mpNumbering(Kratos::make_unique<DistributedNumbering<TIndexType>>(*rOther.mpNumbering)),
        mToRecvLocalId(rOther.mToRecvLocalId),
        mPositionWithinData(rOther.mPositionWithinData),
        mVectorCommColors(rOther.mVectorCommColors)
    {
    }

    ///this function "returns" (writes onto rDestinationVector) a local array containing the values identified by the rGlobalIndices list passed in the constructor
    template< class TDistributedVectorType, class TLocalVectorType, class TApplyFunctorType=std::plus<typename TLocalVectorType::value_type>>
    void Apply(TDistributedVectorType& rDestinationVector, 
               const TLocalVectorType& rLocalDataVector
               ) const
    {
        std::vector<typename TLocalVectorType::value_type> send_buffer;
        std::vector<typename TLocalVectorType::value_type> recv_buffer;
        for(auto color : mVectorCommColors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                const auto& local_ids = mToRecvLocalId.find(color)->second;
                const auto& position_within_data = mPositionWithinData.find(color)->second;
                recv_buffer.resize(local_ids.size());
                send_buffer.resize(0);
                for(IndexType i=0; i<position_within_data.size(); ++i)
                    send_buffer.push_back(rLocalDataVector[position_within_data[i]]);

                //NOTE: this can be made nonblocking
                GetComm().SendRecv(send_buffer, color, 0, recv_buffer, color, 0); 

                //Apply the remotely received data
                for(IndexType i=0; i<recv_buffer.size(); ++i)
                {
                    rDestinationVector[local_ids[i]] = TApplyFunctorType()(rDestinationVector[local_ids[i]],recv_buffer[i]);
                }
            }
        }
        //treat local datas (no communication is needed)
        const auto& local_ids = mToRecvLocalId.find(GetComm().Rank())->second;
        const auto& position_within_data = mPositionWithinData.find(GetComm().Rank())->second;
        for(IndexType i=0; i<position_within_data.size(); ++i)
            rDestinationVector[local_ids[i]] = TApplyFunctorType()( rDestinationVector[local_ids[i]], rLocalDataVector[position_within_data[i]] );
    }

    /// Destructor.
    virtual ~DistributedVectorExporter(){}

    ///@}
    ///@name Operators
    ///@{
    const DataCommunicator& GetComm() const{
        return mrComm;
    }

    
    ///@}
    ///@name Operations
    ///@{

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
    buffer << "DistributedVectorExporter" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DistributedVectorExporter";}

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
    typename DistributedNumbering<IndexType>::UniquePointer mpNumbering;
    std::unordered_map<int, std::vector<IndexType>> mToRecvLocalId;
    std::unordered_map<int, std::vector<IndexType>> mPositionWithinData;
    std::vector<int> mVectorCommColors;

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
    DistributedVectorExporter& operator=(DistributedVectorExporter const& rOther){}

    ///@}

}; // Class DistributedVectorExporter

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                DistributedVectorExporter<TIndexType>& rThis)
                {
                    return rIStream;
                }

/// output stream function
template<class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                const DistributedVectorExporter<TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_VECTOR_EXPORTER_H_INCLUDED  defined


