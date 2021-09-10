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
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_POINTER_MAP_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_POINTER_MAP_COMMUNICATOR_H_INCLUDED


// System includes
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

// External includes
#include "concurrentqueue/concurrentqueue.h"

// Project includes
#include "containers/global_pointers_unordered_map.h"
#include "containers/global_pointers_vector.h"
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/parallel_environment.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
///@addtogroup Kratos MPI Core
///@{

///@name Kratos Classes
///@{

template <class TPointerDataType, class TValueDataType>
class GlobalPointerMapCommunicator;

/**
 * @brief Proxy class to update local and non-local data
 *
 * This class is used to update local and non-local data with given TApplyFunctor and
 * TApplyFunctor is used to update local TPointerDataType data,  Data belonging to remote
 * gps will be stored, and applied after communication in their owning rank.
 *
 * Example:
 * 1. Initial values
 *      Node id   Rank    TEMPERATURE
 *      1         0       10
 *      2         1       20
 *      3         2        5
 *
 * 2. Apply Update method with following Global Pointer map
 *      Rank: 0
 *      Map Key (using node id)         TEMPERATURE
 *      1                               5
 *      2                               10
 *
 *      Rank: 1
 *      Map Key (using node id)         TEMPERATURE
 *      1                               6
 *      3                               2
 *
 *      Rank: 2
 *      Empty Global Pointer map
 *
 * 3. Resulting values after SendAndApplyRemotely method call
 *      Node id   Rank    TEMPERATURE
 *      1         0       10 + 5 + 6 = 21
 *      2         1       20 + 10 = 30
 *      3         2       5 + 2 = 7
 *
 * Signature of TApplyFunctor:
 *      void(TPointerDataType& rPointerDataTypeObject, const TValueDataType& NewValue)
 *
 *
 * @tparam TPointerDataType
 * @tparam TValueDataType
 * @tparam TApplyFunctor
 */
template<class TPointerDataType, class TValueDataType, class TApplyFunctor>
class ApplyProxy
{
public:
    ///@name Type definitions
    ///@{

    using ProxyType = ApplyProxy<TPointerDataType, TValueDataType, TApplyFunctor>;

    using TGPMapCommunicator = GlobalPointerMapCommunicator<TPointerDataType, TValueDataType>;

    using GlobalPointerType = GlobalPointer<TPointerDataType>;

    using DataVectorType = std::vector<TValueDataType>;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Apply Proxy object
     *
     * @param rApplyFunctor                         Thread safe update lambda method
     * @param rPointerCommunicator                  Map communicator
     */
    ApplyProxy(
        const TApplyFunctor& rApplyFunctor,
        TGPMapCommunicator& rPointerCommunicator)
        : mCurrentRank(rPointerCommunicator.GetMyPID()),
          mrApplyFunctor(rApplyFunctor),
          mrPointerCommunicator(rPointerCommunicator)
    {
        KRATOS_TRY

        if (rPointerCommunicator.IsDistributed()) {

            KRATOS_DEBUG_ERROR_IF(OpenMPUtils::IsInParallel() != 0)
                << "Constructing a proxy in a parallel "
                   "region is not allowed.\n";

            this->mUpdateMethod = &ProxyType::AssignLocalAndRemoteData;
        } else {
            this->mUpdateMethod = &ProxyType::AssignLocalData;
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Assigns value of the GlobalPointer
     *
     * This method assigns gp value by rValue using mrApplyFunctor in the case if
     * given rGlobalPointer is a local gp, otherwise mrNonLocalGlobalPointerMapsVector and
     * mrNonLocalDataValueMapsVector are used to store new gp and value,
     * which will be used in mpi communication.
     *
     * In the case of serial job, this method calls AssignLocalData, where no checks are done
     * to ensure gps are local because all gps are local.
     *
     * In case of distributed job, this method calls AssignLocalAndRemoteData where checks are
     * performed to identify gps are local or non-local, based on that mrApplyFunctor
     * methods are called or values corresponding to remote gps are stored.
     *
     * This does not have additional cost in serial run (except for function pointer call)
     *
     * @see GlobalPointerMapCommunicator
     *
     * @param rGP               Global pointer of the destination
     * @param rValue            Value to be used in assigning
     */
    void Assign(
        GlobalPointerType& rGlobalPointer,
        const TValueDataType& rValue)
    {
        (this->*(this->mUpdateMethod))(rGlobalPointer, rValue);
    }

    /**
     * @brief This method does all the communication
     *
     * This method should be called once (or least amount of times) since this involves doing
     * mpi communication to update remote gp values.
     *
     * This does not have any cost in the serial run.
     *
     * Ghost mesh synchronization needs to be done afterwards.
     */
    void SendAndApplyRemotely()
    {
        mrPointerCommunicator.SendAndApplyRemotely(*this);
    }

    ///@}

private:
    ///@name Private members
    ///@{

    const int mCurrentRank;
    const TApplyFunctor& mrApplyFunctor;

    void (ProxyType::*mUpdateMethod)(GlobalPointerType&, const TValueDataType&);

    moodycamel::ConcurrentQueue<std::pair<GlobalPointerType, TValueDataType>> mNonLocalGlobalPointerValueConcurrentQueue;

    TGPMapCommunicator& mrPointerCommunicator;

    ///@}
    ///@name Private operations
    ///@{

    /**
     * @brief Assign local values
     *
     * This method assigns local gp values using mrApplyFunctor.
     * It assumes rGlobalPointer is a local gp always, therefore
     * no checks are performed.
     *
     * @param rGlobalPointer    Local Global pointer of the destination
     * @param rValue            Value to be used in assigning
     */
    void AssignLocalData(
        GlobalPointerType& rGlobalPointer,
        const TValueDataType& rValue)
    {
        mrApplyFunctor(*rGlobalPointer, rValue);
    }


    /**
     * @brief Assign local and non-local values
     *
     * This method assigns local gp values instantly using mrApplyFunctor.
     * As for the non-local gps, they are stored to be used for future
     * communication
     *
     * @param rGlobalPointer    Global pointer of the destination
     * @param rValue            Value to be used in assigning
     */
    void AssignLocalAndRemoteData(
        GlobalPointerType& rGlobalPointer,
        const TValueDataType& rValue)
    {
        const int data_rank = rGlobalPointer.GetRank();

        if (data_rank == mCurrentRank) {
            mrApplyFunctor(*rGlobalPointer, rValue);
        } else {
            mNonLocalGlobalPointerValueConcurrentQueue.enqueue(
                std::make_pair(std::move(rGlobalPointer), rValue));
        }
    }

    ///@}
    ///@name Friend class definitions
    ///@{

    friend class GlobalPointerMapCommunicator<TPointerDataType, TValueDataType>;

    ///@}

};

/// Short class definition.
/** Detail class definition.
*/
template <class TPointerDataType, class TValueDataType>
class GlobalPointerMapCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using GlobalPointerType = GlobalPointer<TPointerDataType>;

    using DataVectorType = std::vector<TValueDataType>;

    using GlobalPointerValueVectorPair = std::pair<GlobalPointerType, DataVectorType>;

    template<class TApplyFunctor>
    using ProxyType = ApplyProxy<TPointerDataType, TValueDataType, TApplyFunctor>;

    /// Pointer definition of GlobalPointerMapCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointerMapCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Global Pointer Map Communicator object
     *
     * This constructor should only be called in non-parallel regions
     *
     * @param rDataCommunicator         Data communicator
     */
    GlobalPointerMapCommunicator(
        const DataCommunicator& rDataCommunicator)
        : mrDataCommunicator(rDataCommunicator),
          mCurrentRank(rDataCommunicator.Rank())
    {
    }

    /// Destructor.
    virtual ~GlobalPointerMapCommunicator() = default;

    /// Assignment constructor
    GlobalPointerMapCommunicator& operator=(GlobalPointerMapCommunicator const& rOther)  = delete;

    /// Copy constructor.
    GlobalPointerMapCommunicator(GlobalPointerMapCommunicator const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Sends and applies remote Gp values
     *
     * This method is used to do communication between processes to apply non_local gp
     * values in each process in GP's owning process
     *
     * Ghost mesh synchronization needs to be done afterwards.
     *
     * @tparam TApplyFunctor
     * @param rApplyProxy           The proxy which holds the user specified methods
     */
    template <class TApplyFunctor>
    void SendAndApplyRemotely(ProxyType<TApplyFunctor>& rApplyProxy)
    {
        KRATOS_TRY

        if (IsDistributed()) {

            KRATOS_DEBUG_ERROR_IF(OpenMPUtils::IsInParallel() != 0)
                << "Calling SendAndApplyRemotely in a parallel region is not "
                   "allowed.\n";

            // get the final map for communications
            // constructing an rank based GlobalPointersUnorderedMap will make unique keys list in GlobalPointersUnorderedMap
            // which will result in lower communication for keys. These unique keys can be easily OMP parallelized
            std::unordered_map<int, GlobalPointersUnorderedMap<TPointerDataType, DataVectorType>> non_local_map;

            bool found_gp_value = true;
            while (found_gp_value) {
                std::pair<GlobalPointerType, TValueDataType> item_pair;
                found_gp_value =
                    rApplyProxy.mNonLocalGlobalPointerValueConcurrentQueue.try_dequeue(item_pair);
                if (found_gp_value) {
                    auto& r_gp = item_pair.first;
                    non_local_map[r_gp.GetRank()][r_gp].push_back(
                        std::move(item_pair.second));
                }
            }

            // compute the communication plan
            const auto& colors = ComputeCommunicationPlan(non_local_map);

            // perform send and receives to get and send remote data
            for (const auto color : colors) {
                if (color >= 0) {
                    // In here
                    //      1. We can do communication twice using default serializer for two map vectors.
                    //              First map: (rank, gp)
                    //              Second map: (rank, data_value)
                    //      2. We can combine non_local_gp_map and non_local_data_map to one map, and do communication once
                    //         after serializing the custom map.
                    //              Current implementation
                    //
                    // Pros and Cons of 1:
                    //      default serialization is used -> Pro
                    //      OMP parallel loops cannot be used because received_gps_vector may not be unique -> Con
                    //          This will be a performance hit if the lambda proxy is computationally expensive
                    //      Communication need to be done twice once for gps, once for values -> Con
                    // Pros and Const of 2:
                    //      custom serialization will be used -> Con
                    //      OMP parallel loops can be used in applying the proxy because the gps will be unique in the map -> Pro
                    //      communication need to be done only once -> Pro
                    //
                    // I think the 2nd options outweighs the first one, so I implemented the second one. Suggestions are
                    // greately appreciated.

                    const auto& received_gp_map = mrDataCommunicator.SendRecv(non_local_map[color], color, color);

                    // create list for OMP parallel looping
                    // using move semantics to move data of GP and std::vector<TValueDataType>
                    // objects in received_gp_map
                    std::vector<GlobalPointerValueVectorPair> gp_value_pair_list;
                    gp_value_pair_list.resize(received_gp_map.size());
                    for (const auto& r_pair : received_gp_map) {
                        gp_value_pair_list.push_back(std::move(r_pair));
                    }

                    // running this in parallel assuming rApplyProxy.mrApplyFunctor is thread safe
                    BlockPartition<std::vector<GlobalPointerValueVectorPair>>(gp_value_pair_list).for_each([&](GlobalPointerValueVectorPair& rItem) {
                        auto& r_pointer_data_type_entity = *(rItem.first);
                        for (IndexType i = 0; i < rItem.second.size(); ++i) {
                            // it is safer to call the serial update method here
                            // because we should only get process's local gps when communication is done
                            rApplyProxy.mrApplyFunctor(r_pointer_data_type_entity, rItem.second[i]);
                        }
                    });
                }
            }

            // Clearing of the concurrent queue is not required here because, when dequeueing, it will make an
            // empty concurrent queue.
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the Apply Proxy object
     *
     * The functor passed via rApplyFunctor should be thread safe
     *
     * Returns the Apply proxy.
     *
     * @tparam TApplyFunctor
     * @param rApplyFunctor                Thread safe functor
     * @return ApplyProxy<TPointerDataType, TValueDataType, TApplyFunctor>
     */
    template <class TApplyFunctor>
    ProxyType<TApplyFunctor> GetApplyProxy(TApplyFunctor&& rApplyFunctor)
    {
        return ProxyType<TApplyFunctor>(std::forward<TApplyFunctor>(rApplyFunctor), *this);
    }

    /**
     * @brief Returns jobs parallel status
     *
     */
    bool IsDistributed() const
    {
        return mrDataCommunicator.IsDistributed();
    }

    /**
     * @brief Get the current rank
     *
     */
    int GetMyPID() const
    {
        return mCurrentRank;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GlobalPointerMapCommunicator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointerMapCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    const DataCommunicator& mrDataCommunicator;
    const int mCurrentRank;

    ///@}
    ///@name Protected Operators
    ///@{

    template<class... TArgs>
    std::vector<int> ComputeCommunicationPlan(const std::unordered_map<int, TArgs...>& rNonLocalGlobalPointerMap)
    {
        std::vector<int> send_list;
        send_list.reserve(rNonLocalGlobalPointerMap.size());
        for (const auto& r_pair : rNonLocalGlobalPointerMap) {
            send_list.push_back(r_pair.first);
        }
        std::sort(send_list.begin(), send_list.end());

        return MPIColoringUtilities::ComputeCommunicationScheduling(
            send_list, mrDataCommunicator);
    }

    ///@}

}; // Class GlobalPointerMapCommunicator

///@}
///@name Input and output
///@{


/// input stream function
template <class TPointerDataType, class TDataType>
inline std::istream& operator>>(
    std::istream& rIStream,
    GlobalPointerMapCommunicator<TPointerDataType, TDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TPointerDataType, class TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const GlobalPointerMapCommunicator<TPointerDataType, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POINTER_MAP_COMMUNICATOR_H_INCLUDED  defined


