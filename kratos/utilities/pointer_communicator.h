//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

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

/**
 * @class ResultsProxy
 * @brief A template class to proxy results, whether they are locally or remotely owned.
 * @details This class acts as a proxy for results. It holds data in a map and applies
 * the user-provided functor to it. If data is locally owned, it applies the functor
 * directly. If the data is remotely owned, it retrieves it from the map.
 * @tparam TPointerDataType  - The datatype for the pointer.
 * @tparam TFunctorType - The functor type to be used in processing.
 */
template< class TPointerDataType, class TFunctorType >
class ResultsProxy
{
public:
    /// Type alias for the result of applying the functor to a global pointer of TPointerDataType
    using TSendType = std::invoke_result_t<TFunctorType,GlobalPointer<TPointerDataType>&>;

    /**
     * @brief Constructor.
     * @param current_rank Current rank of the node.
     * @param NonLocalData Data that is not local to the node.
     * @param UserFunctor Functor to be used for computation.
     * @param pPointerComm Pointer to the communicator.
     */
    ResultsProxy(
        int current_rank,
        GlobalPointersUnorderedMap< TPointerDataType, TSendType > NonLocalData,
        TFunctorType UserFunctor,
        GlobalPointerCommunicator<TPointerDataType>* pPointerComm
    ):
        mCurrentRank(current_rank), mNonLocalData(NonLocalData), mUserFunctor(UserFunctor), mpPointerComm(pPointerComm)
    {}

    /// Destructor.
    virtual ~ResultsProxy() {}

    /**
     * @brief Returns the effect of "user_function(rGlobalPointer)" whether the rGlobalPointer is locally or remotely owned.
     * @param rGlobalPointer The pointer that is being checked.
     * @return The result of applying the functor to the pointer.
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

    /**
     * @brief Overload for constant GlobalPointer.
     * @param rGlobalPointer The constant pointer that is being checked.
     * @return The result of applying the functor to the pointer.
     */
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

    /**
     * @brief Checks if the node has the GlobalPointer.
     * @param rGlobalPointer The pointer that is being checked.
     * @return True if the node has the pointer, otherwise False.
     */
    bool Has(GlobalPointer<TPointerDataType>& rGlobalPointer) const
    {
        if(rGlobalPointer.GetRank() == mCurrentRank)
            return true;
        else
            return mNonLocalData.find(rGlobalPointer) != mNonLocalData.end();
    }

    /**
     * @brief Overload for constant GlobalPointer.
     * @param rGlobalPointer The constant pointer that is being checked.
     * @return True if the node has the pointer, otherwise False.
     */
    bool Has(const GlobalPointer<TPointerDataType>& rGlobalPointer) const
    {
        if(rGlobalPointer.GetRank() == mCurrentRank)
            return true;
        else
            return mNonLocalData.find(rGlobalPointer) != mNonLocalData.end();
    }

    /**
     * @brief Updates the NonLocalData using the user functor and communicator.
     */
    void Update()
    {
        mpPointerComm->Update(mUserFunctor, mNonLocalData);
    }

private:
    const int mCurrentRank;  // The current rank of the node.
    GlobalPointersUnorderedMap< TPointerDataType, TSendType > mNonLocalData;  // Map holding data not local to the node.
    TFunctorType mUserFunctor;  // Functor to be used for computation.
    GlobalPointerCommunicator<TPointerDataType>* mpPointerComm;  // Pointer to the communicator.

};


/**
 * @class GlobalPointerCommunicator
 * @brief A template class for handling communication related to global pointers.
 * @details This class manages communication related to global pointers of a specified type.
 * It is responsible for creating a communication plan for sending/receiving global pointers,
 * and for applying a user-provided function to global pointers.
 * @tparam TPointerDataType The datatype for the pointer.
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

    /**
     * @brief Default constructor.
     * @param rComm The data communicator.
     * @param rGpList List of global pointers to be added.
     */
    GlobalPointerCommunicator(const DataCommunicator& rComm, GlobalPointersVector< TPointerDataType >& rGpList ):
        mrDataCommunicator(rComm)
    {
        AddPointers(rGpList.ptr_begin(),rGpList.ptr_end());

        if(mrDataCommunicator.IsDistributed())
        {
            ComputeCommunicationPlan();
        }
    }

    /**
     * @brief Constructor with iterator range for global pointers.
     * @tparam TIteratorType Iterator type.
     * @param rComm The data communicator.
     * @param itBegin Beginning of the range.
     * @param itEnd End of the range.
     */
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

    /**
     * @brief Constructor with functor for generating global pointers.
     * @tparam TFunctorType Functor type.
     * @param rComm The data communicator.
     * @param rFunctor The functor for generating global pointers.
     */
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

    /**
     * @brief Constructor with functor and default data communicator.
     * @tparam TFunctorType Functor type.
     * @param rFunctor The functor for generating global pointers.
     */
    template< class TFunctorType >
    GlobalPointerCommunicator(TFunctorType rFunctor)
        : GlobalPointerCommunicator(ParallelEnvironment::GetDefaultDataCommunicator(), rFunctor)
    {}

    /// Destructor.
    virtual ~GlobalPointerCommunicator() {}

    /**
     * @brief Applies a user-provided function to the global pointers and return a proxy to the results.
     * @tparam TFunctorType Functor type.
     * @param UserFunctor The user-provided function.
     * @return A proxy to the results.
     */
    template< class TFunctorType >
    ResultsProxy<
    TPointerDataType,
    TFunctorType // TODO: Unfortunately this is deprecated in c++17, so we will have to change this call in the future
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

    /**
     * @brief Updates the non-local data using a user-provided function.
     * @tparam TFunctorType Functor type.
     * @param rUserFunctor The user-provided function.
     * @param rNonLocalData Non-local data to be updated.
     */
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

    /**
     * @brief Returns the data communicator.
     * @return The data communicator.
     */
    const DataCommunicator& GetDataCommunicator() const
    {
        return mrDataCommunicator;
    }

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

    std::unordered_map<int, GlobalPointersVector< TPointerDataType > > mNonLocalPointers; /// Non local pointers stored as map with rank as key.

    const DataCommunicator& mrDataCommunicator; /// DataCommunicator reference.

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * @brief Adds global pointers to the communicator.
     * @details This function adds the global pointers from the provided range to the communicator.
     * It only adds the pointers that are not on the current rank. After adding the pointers,
     * it ensures that no repeated information is stored.
     * @tparam TIteratorType Iterator type.
     * @param begin Beginning of the range of global pointers.
     * @param end End of the range of global pointers.
     */
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

    /**
     * @brief Computes the communication plan for the communicator.
     * @details This function computes the communication plan based on the ranks of the non-local pointers.
     * It determines the scheduling for the communication, which is stored in the 'mColors' member variable.
     */
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

    std::vector<int> mColors; /// Communication colors (indicating ranks involved in communication)

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


