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

/// Short class definition.
/** Detail class definition.
*/
template< class TPointerDataType, class TSendType >
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
    GlobalPointerCommunicator(DataCommunicator& rComm)
    mrDataCommunicator(rComm)
    {}

    /// Destructor.
    virtual ~GlobalPointerCommunicator(){}

    void AddPointer( GlobalPointer<TPointerDataType>& gp)
    {
        if(gp.GetRank() != mrDataCommunicator.Rank())
            mNonLocalPointers[p->GetRank()].push_back(gp);
    }

    /**this function applies the input function onto the remote data
     * and retrieves the result to the caller rank
     */
    void Apply(std::function< TSendType(GlobalPointer<TPointerDataType>&) >& function)
    {
        //TODO: avoid doing anything if not distributed
        
        mfunction = function; //storing the function for later use
        mNonLocalData.clear();

        //compute communication plan
        std::vector<int> send_list.reserve( mNonLocalPointers.size() );
        for(auto& it : NonLocalPointers)
            send_list.push_back( it->first );
        std::sort(send_list.begin(), send_list.end())
        auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list[current_rank], mrDataCommunicator);
        int send_rank = mrDataCommunicator.Rank();

        //sendrecv data
        for(auto color : colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                auto& gps_to_be_sent = mNonLocalPointers[color];

                //TODO: pass Unique to the mNonLocalPointers[recv_rank]

                auto recv_global_pointers = GenericSendRecv(gps_to_be_sent, color, send_rank );

                std::vector< TSendType > locally_gathered_data; //this is local but needs to be sent to the remote node
                for(auto& gp : recv_global_pointers)
                    locally_gathered_data.push_back( mfunction(gp) );

                auto remote_data = GenericSendRecv(locally_gathered_data, color, send_rank );

                for(unsigned int i=0; i<remote_data.size(); ++i)
                    mNonLocalData[gps_to_be_sent[i]] = remote_data[i];
            }
        }
    }


    /**this function returns the effect of "function(gp)" both if the gp is locally owned 
     * and if it is remotely owned
     */    
    TSendType Get(GlobalPointer<TPointerDataType>& gp);
    {
        if(gp.GetRank() == mrDataCommunicator.Rank())
        {
            return mfunction(gp);
        }
        else
        {
            return mNonLocalData(gp);
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
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GlobalPointerCommunicator";}

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
    GlobalPointersMap< TPointerDataType, TSendType > mNonLocalData;
    DataCommunicator mrDataCommunicator;
    std::function< TSendType(GlobalPointer<TPointerDataType>&) > mfunction;

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


    ///@}
    ///@name Private Operators
    ///@{
    template< TDataType>
    TDataType GenericSendRecv(TDataType& send_buffer, int send_rank, int recv_rank)
    {
        MpiSerializer send_serializer;
        recv_serializer.save("data",send_buffer);
        std::string send_buffer = send_serializer.GetStringRepresentation();

        std::string recv_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank);

        MpiSerializer recv_serializer;

        TDataType recv_data(recv_buffer);
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
    GlobalPointerCommunicator& operator=(GlobalPointerCommunicator const& rOther){}

    /// Copy constructor.
    GlobalPointerCommunicator(GlobalPointerCommunicator const& rOther){}

    ///@}

}; // Class GlobalPointerCommunicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                GlobalPointerCommunicator& rThis){}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const GlobalPointerCommunicator& rThis)
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


