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


#if !defined(KRATOS_POINTER_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_POINTER_COMMUNICATOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"
#include "includes/global_pointer.h"
#include "includes/mpi_serializer.h"
#include "containers/global_pointers_vector.h"
#include "containers/global_pointers_unordered_map.h"

#include "utilities/mpi_coloring_utilities.h"

namespace Kratos
{
///@addtogroup Kratos MPI Core
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
template< class TPointerDataType, class TSendType >
class ResultsProxy
{
public:

    ResultsProxy(
        int current_rank,
        GlobalPointersUnorderedMap< TPointerDataType, TSendType > NonLocalData,
        std::function< TSendType(GlobalPointer<TPointerDataType>&) > user_function
    ):
        mCurrentRank(current_rank), mNonLocalData(NonLocalData), mUserFunction(user_function)
    {}

    /// Destructor.
    virtual ~ResultsProxy() {}

    /**this function returns the effect of "user_function(gp)" both if the gp is locally owned
     * and if it is remotely owned
     */
    TSendType Get(GlobalPointer<TPointerDataType>& gp)
    {
        if(gp.GetRank() == mCurrentRank)
            return mUserFunction(gp);
        else
            return mNonLocalData[gp];
    }

private:
    const int mCurrentRank;
    GlobalPointersUnorderedMap< TPointerDataType, TSendType > mNonLocalData;
    std::function< TSendType(GlobalPointer<TPointerDataType>&) > mUserFunction;

};


/// Short class definition.
/** Detail class definition.
*/
template< class TPointerDataType >
class GlobalPointerCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GlobalPointerCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointerCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GlobalPointerCommunicator(const DataCommunicator& rComm, GlobalPointersVector< TPointerDataType >& gp_list ):
        mrDataCommunicator(rComm)
    {
        AddPointers(gp_list.begin(),gp_list.end());

        if(mrDataCommunicator.IsDistributed())
        {
            //compute communication plan
            std::vector<int> send_list;
            send_list.reserve( mNonLocalPointers.size() );
            for(auto& it : mNonLocalPointers)
                send_list.push_back( it.first );

            std::sort(send_list.begin(), send_list.end());
            mColors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, mrDataCommunicator);
        }
    }

    template< class TIteratorType >
    GlobalPointerCommunicator(const DataCommunicator& rComm, TIteratorType begin, TIteratorType end):
        mrDataCommunicator(rComm)
    {
        AddPointers(begin,end);

        if(mrDataCommunicator.IsDistributed())
        {
            //compute communication plan
            std::vector<int> send_list;
            send_list.reserve( mNonLocalPointers.size() );
            for(auto& it : mNonLocalPointers)
                send_list.push_back( it.first );

            std::sort(send_list.begin(), send_list.end());
            mColors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, mrDataCommunicator);
        }
    }

    /// Destructor.
    virtual ~GlobalPointerCommunicator() {}

    /**this function applies the input function onto the remote data
     * and retrieves the result to the caller rank
     */
    template< class TSendType >
    ResultsProxy<TPointerDataType, TSendType> Apply(std::function< TSendType(GlobalPointer<TPointerDataType>&) > user_function)
    {
        if(mrDataCommunicator.IsDistributed())
        {
            const int current_rank = mrDataCommunicator.Rank();

            GlobalPointersUnorderedMap< TPointerDataType, TSendType > non_local_data;


            //sendrecv data
            for(auto color : mColors)
            {
                if(color >= 0) //-1 would imply no communication
                {
                    auto& gps_to_be_sent = mNonLocalPointers[color];

                    //TODO: pass Unique to the mNonLocalPointers[recv_rank]

                    auto recv_global_pointers = SendRecv(gps_to_be_sent, color, color );

                    std::vector< TSendType > locally_gathered_data; //this is local but needs to be sent to the remote node
                    for(auto& gp : recv_global_pointers)
                        locally_gathered_data.push_back( user_function(gp) );

                    auto remote_data = SendRecv(locally_gathered_data, color, color );

                    for(unsigned int i=0; i<remote_data.size(); ++i)
                        non_local_data[gps_to_be_sent[i]] = remote_data[i];
                }
            }

            return ResultsProxy<TPointerDataType, TSendType>(current_rank,non_local_data,user_function );
        }
        else
        {
            GlobalPointersUnorderedMap< TPointerDataType, TSendType > non_local_data;
            return ResultsProxy<TPointerDataType, TSendType>(0,non_local_data,user_function );
        }

    }
    ///@}
    ///@name Operators
    ///@{


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
        buffer << "GlobalPointerCommunicator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointerCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    std::unordered_map<int, GlobalPointersVector< TPointerDataType > > mNonLocalPointers;
    const DataCommunicator& mrDataCommunicator;


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{
    template< class TIteratorType >
    void AddPointers( TIteratorType begin, TIteratorType end)
    {
        if(mrDataCommunicator.IsDistributed())
        {
            const int current_rank = mrDataCommunicator.Rank();
            for(auto it = begin; it != end; ++it)
                if(it->GetRank() != current_rank)
                    mNonLocalPointers[it->GetRank()].push_back(*it);
        }
    }

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
    std::vector<int> mColors;


    ///@}
    ///@name Private Operators
    ///@{
    template< class TDataType>
    TDataType SendRecv(TDataType& send_buffer, int send_rank, int recv_rank)
    {
        MpiSerializer send_serializer;
        send_serializer.save("data",send_buffer);
        std::string send_string = send_serializer.GetStringRepresentation();

        std::string recv_string = mrDataCommunicator.SendRecv(send_string, send_rank, send_rank);

        MpiSerializer recv_serializer(recv_string);

        TDataType recv_data;
        recv_serializer.load("data",recv_data);
        return recv_data;
    }

    //TODO explicitly instantiate this function to SendRecv for basic types


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
    GlobalPointerCommunicator& operator=(GlobalPointerCommunicator const& rOther) {}

    /// Copy constructor.
    GlobalPointerCommunicator(GlobalPointerCommunicator const& rOther) {}

    ///@}

}; // Class GlobalPointerCommunicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TPointerDataType, class TSendType >
inline std::istream& operator >> (std::istream& rIStream,
                                  GlobalPointerCommunicator<TPointerDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TPointerDataType, class TSendType >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GlobalPointerCommunicator<TPointerDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POINTER_COMMUNICATOR_H_INCLUDED  defined


