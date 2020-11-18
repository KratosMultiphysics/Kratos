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

// External includes

// Project includes
#include "containers/global_pointers_unordered_map.h"
#include "containers/global_pointers_vector.h"
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/parallel_environment.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/parallel_utilities.h"

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
 * This class is used to update local and non-local data with given TLocalApplyFunctor and TNonLocalApplyFunctor
 * TLocalApplyFunctor is used to update local TPointerDataType data, TNonLocalApplyFunctor is used to update
 * non local global pointer map data which then will use TLocalApplyFunctor after communication to update value
 * at the owner rank.
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
 * Signature of TLocalApplyFunctor:
 *      void(TPointerDataType& rPointerDataTypeObject, const TValueDataType& NewValue)
 *
 * Signature of TNonLocalApplyFunctor:
 *      void(TValueDataType& rCurrentValueInNonLocalGlobalPointerMap, const TValueDataType& NewValue)
 *
 * @tparam TPointerDataType
 * @tparam TValueDataType
 * @tparam TLocalApplyFunctor
 * @tparam TNonLocalApplyFunctor
 */
template<class TPointerDataType, class TValueDataType, class TLocalApplyFunctor, class TNonLocalApplyFunctor>
class ApplyProxy
{
public:
    ///@name Type definitions
    ///@{

    using ProxyType = ApplyProxy<TPointerDataType, TValueDataType, TLocalApplyFunctor, TNonLocalApplyFunctor>;

    using TGPVector = GlobalPointersVector<TPointerDataType>;

    using TGPDataMap = GlobalPointersUnorderedMap<TPointerDataType, TValueDataType>;

    using TGPNonLocalDataMap = std::unordered_map<int, TGPDataMap>;

    using TGPMapCommunicator = GlobalPointerMapCommunicator<TPointerDataType, TValueDataType>;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Apply Proxy object
     *
     * @param CurrentRank               Current rank
     * @param rLocalApplyFunctor        Thread safe update lambda method
     * @param rNonLocalApplyFunctor     Thread safe non-local gp map value update method
     * @param rNonLocalDataMap          Non-local gp map
     * @param rPointerCommunicator      Map communicator
     */
    ApplyProxy(
        const int CurrentRank,
        const TLocalApplyFunctor& rLocalApplyFunctor,
        const TNonLocalApplyFunctor& rNonLocalApplyFunctor,
        TGPNonLocalDataMap& rNonLocalDataMap,
        TGPMapCommunicator& rPointerCommunicator)
        : mCurrentRank(CurrentRank),
          mrLocalApplyFunctor(rLocalApplyFunctor),
          mrNonLocalApplyFunctor(rNonLocalApplyFunctor),
          mrNonLocalDataMap(rNonLocalDataMap),
          mrPointerCommunicator(rPointerCommunicator)
    {
        // identify proper methods to be used in updates
        if (mrPointerCommunicator.IsDistributed()) {
            this->mUpdateMethod = &ProxyType::UpdateLocalAndRemoteData;
        } else {
            this->mUpdateMethod = &ProxyType::UpdateLocalData;
        }
    }

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Updates values with given values map
     *
     * This method updates local gps using mrLocalApplyFunctor instantly. Value updates
     * belonging to remote gps are stored in the mrNonLocalDataMap using mrNonLocalApplyFunctor.
     *
     * In the case of serial job, this method calls UpdateLocalData, where no checks are done
     * to ensure gps are local because all gps are local.
     *
     * In case of distributed job, this method calls UpdateLocalAndRemoteData where checks are
     * performed to identify gps are local or non-local, based on that mrLocalApplyFunctor or mrNonLocalApplyFunctor
     * methods are called.
     *
     * The input rGPDataMap map should not contain any key gps which are not included in creating the
     * GlobalPointerMapCommunicator used in this proxy (i.e. mrPointerCommunicator)
     *
     * This method can be called several times without severe additional computational cost
     * since this does not do any mpi communication.
     *
     * This does not have additional cost in serial run
     *
     * @see GlobalPointerMapCommunicator
     *
     * @param rGPDataMap        Input values map (key: GlobalPointer<TPointerDataType>, value: TValueDataType)
     */
    void Update(const TGPDataMap& rGPDataMap)
    {
        // get gp vector
        const auto& gps = TGPMapCommunicator::GetKeys(rGPDataMap);

        // running this in parallel assuming mrLocalApplyFunctor is thread safe
        IndexPartition<int>(gps.size()).for_each([&](const int Index) {
            auto p_itr = rGPDataMap.find(gps[Index]);
            Update(p_itr->first, p_itr->second);
        });
    }

    /**
     * @brief Update value of the GlobalPointer
     *
     * This method updates gp value by rNewValue using mrLocalApplyFunctor in the case if
     * given rGlobalPointer is a local gp, otherwise mrNonLocalDataMap is updated using
     * mrNonLocalApplyFunctor.
     *
     * @param rGP               Global pointer of the destination
     * @param rNewValue         New value to be used in update
     */
    void Update(
        const GlobalPointer<TPointerDataType>& rGlobalPointer,
        const TValueDataType& rNewValue)
    {
        (this->*(this->mUpdateMethod))(rGlobalPointer, rNewValue);
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
     *
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
    const TLocalApplyFunctor& mrLocalApplyFunctor;
    const TNonLocalApplyFunctor& mrNonLocalApplyFunctor;

    void (ProxyType::*mUpdateMethod)(const GlobalPointer<TPointerDataType>&, const TValueDataType&);

    TGPNonLocalDataMap& mrNonLocalDataMap;
    TGPMapCommunicator& mrPointerCommunicator;

    ///@}
    ///@name Private operations
    ///@{

    /**
     * @brief Update local values
     *
     * This method updates local gp values using mrLocalApplyFunctor.
     * It assumes rGlobalPointer is a local gp always, therefore
     * no checks are performed.
     *
     * @param rGlobalPointer    Local Global pointer of the destination
     * @param rNewValue         New value to be used in update
     */
    void UpdateLocalData(
        const GlobalPointer<TPointerDataType>& rGlobalPointer,
        const TValueDataType& rNewValue)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rGlobalPointer.GetRank() != mCurrentRank)
            << "Using local global pointer update method with a non-local "
               "global pointer. [ MyPID = "
            << mCurrentRank
            << " GlobalPointerRank = " << rGlobalPointer.GetRank() << " ].\n";

        auto gp_pointer = rGlobalPointer;  // required to remove the const from rGlobalPointer
        mrLocalApplyFunctor(*gp_pointer, rNewValue);

        KRATOS_CATCH("");
    }


    /**
     * @brief Update local and non-local values
     *
     * This method updates local gp values instantly using mrLocalApplyFunctor.
     * As for the non-local gps, it uses mrNonLocalApplyFunctor to update
     * mrNonLocalDataMap values which is used for communication
     *
     * @param rGlobalPointer    Global pointer of the destination
     * @param rNewValue         New value to be used in update
     */
    void UpdateLocalAndRemoteData(
        const GlobalPointer<TPointerDataType>& rGlobalPointer,
        const TValueDataType& rNewValue)
    {
        KRATOS_TRY

        const int data_rank = rGlobalPointer.GetRank();

        if (data_rank == mCurrentRank) {
            UpdateLocalData(rGlobalPointer, rNewValue);
        } else {
            // update mrNonLocalDataMap map for remote data,
            // which will be used in future for communication
            auto p_non_local_map = mrNonLocalDataMap.find(data_rank);

            KRATOS_DEBUG_ERROR_IF(p_non_local_map == mrNonLocalDataMap.end())
                << "Data rank not found in non local maps. [ MyPID = " << mCurrentRank
                << ", DataRank = " << data_rank << " ].\n";

            auto p_non_local_itr = p_non_local_map->second.find(rGlobalPointer);

            KRATOS_DEBUG_ERROR_IF(p_non_local_itr == p_non_local_map->second.end())
                << "Global pointer not found in rank " << data_rank
                << " non local map. [ MyPID = " << mCurrentRank << " ].\n";

            mrNonLocalApplyFunctor(p_non_local_itr->second, rNewValue);
        }

        KRATOS_CATCH("");
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

    template<class TLocalApplyFunctor, class TNonLocalApplyFunctor>
    using ProxyType = ApplyProxy<TPointerDataType, TValueDataType, TLocalApplyFunctor, TNonLocalApplyFunctor>;

    using TGPVector = GlobalPointersVector<TPointerDataType>;

    using TGPDataMap = GlobalPointersUnorderedMap<TPointerDataType, TValueDataType>;

    using TGPNonLocalDataMap = std::unordered_map<int, TGPDataMap>;

    /// Pointer definition of GlobalPointerMapCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointerMapCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Global Pointer Map Communicator object
     *
     * This constructor can be used if there is already a GPVector
     * constructed by the user.
     *
     * In parallel, this constructor will compute communication scheduling
     * for given list of GPs.
     *
     * The rInitializationValue is used to initialize the gp map created
     * for GPs provided by rInitializationValue
     *
     * @param rDataCommunicator         Data communicator
     * @param rGPVector                 Global pointers vector
     * @param rInitializationValue      Initialization value
     */
    GlobalPointerMapCommunicator(
        const DataCommunicator& rDataCommunicator,
        const TGPVector& rGPVector,
        const TValueDataType& rInitializationValue)
        : mrDataCommunicator(rDataCommunicator),
          mCurrentRank(rDataCommunicator.Rank())
    {
        if (IsDistributed()) {
            AddPointers(rGPVector, rInitializationValue);
            ComputeCommunicationPlan();
        }
    }

    /**
     * @brief Construct a new Global Pointer Map Communicator object
     *
     * This constructor uses a functor to create the Gps vector
     *
     * @tparam TVectorFunctorType
     * @param rDataCommunicator         Data communicator
     * @param rInitializationValue      Initialization value
     * @param rVectorFunctor            Functor which returns a GPVector
     */
    template <class TVectorFunctorType>
    GlobalPointerMapCommunicator(
        const DataCommunicator& rDataCommunicator,
        const TValueDataType& rInitializationValue,
        TVectorFunctorType&& rVectorFunctor)
        : mrDataCommunicator(rDataCommunicator),
          mCurrentRank(rDataCommunicator.Rank())
    {
        if (IsDistributed()) {
            const auto& r_gp_vector = rVectorFunctor(mrDataCommunicator);
            AddPointers(r_gp_vector, rInitializationValue);
            ComputeCommunicationPlan();
        }
    }

    template <class TVectorFunctorType>
    GlobalPointerMapCommunicator(
        const TValueDataType& rInitializationValue,
        TVectorFunctorType&& rFunctor)
        : GlobalPointerMapCommunicator(
              ParallelEnvironment::GetDefaultDataCommunicator(),
              rInitializationValue,
              std::forward<TVectorFunctorType>(rFunctor))
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
     * @tparam TLocalApplyFunctor
     * @tparam TNonLocalApplyFunctor
     * @param rApplyProxy           The proxy which holds the user specified methods
     */
    template <class TLocalApplyFunctor, class TNonLocalApplyFunctor>
    void SendAndApplyRemotely(
        ProxyType<TLocalApplyFunctor, TNonLocalApplyFunctor>& rApplyProxy)
    {
        if (IsDistributed()) {
            for (auto color : mColors) {
                if (color >= 0) {
                    const auto& received_gp_map = mrDataCommunicator.SendRecv(
                        mrNonLocalPointers[color], color, color);

                    // get gp vector
                    const auto& gps = GetKeys(received_gp_map);

                    // running this in parallel assuming mrLocalApplyFunctor is thread safe
                    IndexPartition<int>(gps.size()).for_each([&](const int Index) {
                        auto p_itr = received_gp_map.find(gps[Index]);
                        // it is safer to call the serial update method here
                        // because we should only get process's local gps when communication is done
                        rApplyProxy.UpdateLocalData(p_itr->first, p_itr->second);
                    });
                }
            }
        }
    }

    /**
     * @brief Get the Apply Proxy object
     *
     * Returns the Apply proxy.
     *
     * @tparam TLocalApplyFunctor
     * @tparam TNonLocalApplyFunctor
     * @param rApplyFunctor
     * @return ApplyProxy<TPointerDataType, TValueDataType, TApplyFunctor>
     */
    template <class TLocalApplyFunctor, class TNonLocalApplyFunctor>
    ProxyType<TLocalApplyFunctor, TNonLocalApplyFunctor> GetApplyProxy(
        TLocalApplyFunctor&& rLocalApplyFunctor,
        TNonLocalApplyFunctor&& rNonLocalApplyFunctor)
    {
        return ProxyType<TLocalApplyFunctor, TNonLocalApplyFunctor>(
            mCurrentRank, std::forward<TLocalApplyFunctor>(rLocalApplyFunctor),
            std::forward<TNonLocalApplyFunctor>(rNonLocalApplyFunctor),
            mrNonLocalPointers, *this);
    }

    /**
     * @brief Returns jobs parallel status
     *
     */
    bool IsDistributed() const
    {
        return mrDataCommunicator.IsDistributed();
    }

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Get keys vector from a map
     *
     * @tparam TKey                 Key type
     * @tparam TArgs                Additional args
     * @param rMap                  Map input
     * @return std::vector<TKey>    Vector of keys
     */
    template<class TKey, class... TArgs>
    static std::vector<TKey> GetKeys(const std::unordered_map<TKey, TArgs...>& rMap)
    {
        std::vector<TKey> keys;
        keys.resize(rMap.size());

        int local_index = 0;
        for (const auto& r_item : rMap) {
            keys[local_index++] = r_item.first;
        }

        return keys;
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

    TGPNonLocalDataMap mrNonLocalPointers;
    const DataCommunicator& mrDataCommunicator;
    const int mCurrentRank;

    ///@}
    ///@name Protected Operators
    ///@{

    void AddPointers(
        const TGPVector& rGPVector,
        const TValueDataType& rInitializationValue)
    {
        // this method is called always at the constructor, and it is not allowed to have PointerMapCommunicators
        // created in OMP regions so, it is safer to do OMP parallelization here.
        // this should only be called in MPI
        IndexPartition<int>(rGPVector.size()).for_each([&](const int Index) {
            const auto& r_gp = rGPVector(Index);
            if (r_gp.GetRank() != mCurrentRank) {
#pragma omp critical
                {
                    auto& rank_local_map = mrNonLocalPointers[r_gp.GetRank()];
                    rank_local_map[r_gp] = rInitializationValue;
                }
            }
        });
    }

    void ComputeCommunicationPlan()
    {
        auto send_list = GetKeys(mrNonLocalPointers);
        std::sort(send_list.begin(), send_list.end());
        mColors = MPIColoringUtilities::ComputeCommunicationScheduling(
            send_list, mrDataCommunicator);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<int> mColors;

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


