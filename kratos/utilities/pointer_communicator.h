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
#include <type_traits>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"
#include "includes/global_pointer.h"
#include "includes/mpi_serializer.h"
#include "containers/global_pointers_vector.h"
#include "containers/global_pointers_unordered_map.h"
#include "includes/parallel_environment.h"
#include "utilities/communication_coloring_utilities.h"

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
template< class TPointerDataType >
class GlobalPointerCommunicator; //fully defined at the end of the file

template< class TPointerDataType, class TFunctorType >
class ResultsProxy
{
public:

    typedef typename std::result_of< TFunctorType(GlobalPointer<TPointerDataType>&)>::type TSendType;

    ResultsProxy(
        int current_rank,
        GlobalPointersUnorderedMap< TPointerDataType, TSendType > NonLocalData,
        TFunctorType user_functor,
        GlobalPointerCommunicator<TPointerDataType>* pPointerComm
    ):
        mCurrentRank(current_rank), mNonLocalData(NonLocalData), mUserFunctor(user_functor), mpPointerComm(pPointerComm)
    {}

    /// Destructor.
    virtual ~ResultsProxy() {}

    /**this function returns the effect of "user_function(rGlobalPointer)" both if the rGlobalPointer is locally owned
     * and if it is remotely owned
     */
    TSendType Get(GlobalPointer<TPointerDataType>& rGlobalPointer) const
    {
        if(rGlobalPointer.GetRank() == mCurrentRank)
            return mUserFunctor(rGlobalPointer);
        else {
            auto non_local_gp = mNonLocalData.find(rGlobalPointer);
            KRATOS_DEBUG_ERROR_IF(non_local_gp == mNonLocalData.end()) << "Missing entry in NonLocalData" << std::endl;
            return non_local_gp->second;
        }
    }

    TSendType Get(const GlobalPointer<TPointerDataType>& rGlobalPointer) const
    {
        if(rGlobalPointer.GetRank() == mCurrentRank)
            return mUserFunctor(rGlobalPointer);
        else {
            auto non_local_gp = mNonLocalData.find(rGlobalPointer);
            KRATOS_DEBUG_ERROR_IF(non_local_gp == mNonLocalData.end()) << "Missing entry in NonLocalData" << std::endl;
            return non_local_gp->second;
        }
    }
    
    bool Has(GlobalPointer<TPointerDataType>& rGlobalPointer) const 
    {
        if(rGlobalPointer.GetRank() == mCurrentRank)
            return true;
        else
            return mNonLocalData.find(rGlobalPointer) != mNonLocalData.end();
    }   

    bool Has(const GlobalPointer<TPointerDataType>& rGlobalPointer) const
    {
        if(rGlobalPointer.GetRank() == mCurrentRank)
            return true;
        else
            return mNonLocalData.find(rGlobalPointer) != mNonLocalData.end();
    }

    void Update()
    {
        mpPointerComm->Update(mUserFunctor, mNonLocalData);
    }

private:
    const int mCurrentRank;
    GlobalPointersUnorderedMap< TPointerDataType, TSendType > mNonLocalData;
    TFunctorType mUserFunctor;
    GlobalPointerCommunicator<TPointerDataType>* mpPointerComm;

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
    GlobalPointerCommunicator(const DataCommunicator& rComm, GlobalPointersVector< TPointerDataType >& rGpList ):
        mrDataCommunicator(rComm)
    {
        AddPointers(rGpList.ptr_begin(),rGpList.ptr_end());

        if(mrDataCommunicator.IsDistributed())
        {
            ComputeCommunicationPlan();
        }
    }

    template< class TIteratorType >
    GlobalPointerCommunicator(const DataCommunicator& rComm, TIteratorType itBegin, TIteratorType itEnd):
        mrDataCommunicator(rComm)
    {
        AddPointers(itBegin,itEnd);

        if(mrDataCommunicator.IsDistributed())
        {
            ComputeCommunicationPlan();
        }
    }

    template< class TFunctorType >
    GlobalPointerCommunicator(const DataCommunicator& rComm, TFunctorType rFunctor):
        mrDataCommunicator(rComm)
    {
        if(rComm.IsDistributed())
        {
            auto gps = rFunctor(rComm);
            AddPointers(gps.ptr_begin(), gps.ptr_end());
            ComputeCommunicationPlan();
        }
    }

    template< class TFunctorType >
    GlobalPointerCommunicator(TFunctorType rFunctor)
        : GlobalPointerCommunicator(ParallelEnvironment::GetDefaultDataCommunicator(), rFunctor)
    {}

    /// Destructor.
    virtual ~GlobalPointerCommunicator() {}

    /**this function applies the input function onto the remote data
     * and retrieves the result to the caller rank
     */
    template< class TFunctorType >
    ResultsProxy<
    TPointerDataType,
    TFunctorType //unfortunately this is deprecated in c++17, so we will have to change this call in the future
    > Apply(TFunctorType&& UserFunctor)
    {
        typedef typename ResultsProxy<TPointerDataType, TFunctorType >::TSendType SendType;

        const int current_rank = mrDataCommunicator.Rank();

        GlobalPointersUnorderedMap< TPointerDataType, SendType > non_local_data;

        if(mrDataCommunicator.IsDistributed())
        {
            Update(UserFunctor, non_local_data);
        }

        return ResultsProxy<TPointerDataType, TFunctorType>(current_rank,non_local_data,UserFunctor, this );
    }

    template< class TFunctorType >
    void Update(
        TFunctorType& rUserFunctor,
        GlobalPointersUnorderedMap< TPointerDataType, typename ResultsProxy<TPointerDataType, TFunctorType >::TSendType >& rNonLocalData)
    {
        //sendrecv data
        for(auto color : mColors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                auto& gps_to_be_sent = mNonLocalPointers[color];

                auto recv_global_pointers = mrDataCommunicator.SendRecv(gps_to_be_sent, color, color );

                std::vector< typename ResultsProxy<TPointerDataType, TFunctorType >::TSendType > locally_gathered_data; //this is local but needs to be sent to the remote node
                for(auto& gp : recv_global_pointers.GetContainer())
                    locally_gathered_data.push_back( rUserFunctor(gp) );

                auto remote_data = mrDataCommunicator.SendRecv(locally_gathered_data, color, color );

                for(unsigned int i=0; i<remote_data.size(); ++i)
                    rNonLocalData[gps_to_be_sent(i)] = remote_data[i];
            }
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
            {
                auto& gp = *it;
                if(gp.GetRank() != current_rank)
                {
                    mNonLocalPointers[gp.GetRank()].push_back(gp);
                }
            }

            //ensure that no repeated info is stored
            for(auto& non_local : mNonLocalPointers)
                non_local.second.Unique();
        }
    }

    void ComputeCommunicationPlan()
    {
        std::vector<int> send_list;
        send_list.reserve( mNonLocalPointers.size() );
        for(auto& it : mNonLocalPointers)
            send_list.push_back( it.first );

        std::sort(send_list.begin(), send_list.end());
        mColors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, mrDataCommunicator);
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


