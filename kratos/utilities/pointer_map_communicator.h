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
     * @param rApplyFunctor             Thread safe update lambda method
     * @param rNonLocalDataMap          Non-local gp map
     * @param rPointerCommunicator      Map communicator
     */
    ApplyProxy(
        const TApplyFunctor& rApplyFunctor,
        TGPNonLocalDataMap& rNonLocalDataMap,
        TGPMapCommunicator& rPointerCommunicator)
        : mCurrentRank(rPointerCommunicator.GetMyPID()),
          mrApplyFunctor(rApplyFunctor),
          mrNonLocalDataMap(rNonLocalDataMap),
          mrPointerCommunicator(rPointerCommunicator)
    {
        // identify proper methods to be used in updates
        if (mrPointerCommunicator.IsDistributed()) {
            this->mUpdateMethod = &ProxyType::AssignLocalAndRemoteData;
        } else {
            this->mUpdateMethod = &ProxyType::AssignLocalData;
        }

        // identify proper methods to be used in updating non-local gp data
        if (mrPointerCommunicator.IsAssembly()) {
            this->mNonLocalGPDataUpdateMethod = &ProxyType::NonLocalDataAssemblyMethod;
        } else {
            this->mNonLocalGPDataUpdateMethod = &ProxyType::NonLocalDataAssignMethod;
        }
    }

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Assigns values with given values map
     *
     * This method assigns local gps using mrApplyFunctor instantly. Values
     * belonging to remote gps are stored in the mrNonLocalDataMap.
     *
     * In the case of serial job, this method calls AssignLocalData, where no checks are done
     * to ensure gps are local because all gps are local.
     *
     * In case of distributed job, this method calls AssignLocalAndRemoteData where checks are
     * performed to identify gps are local or non-local, based on that mrApplyFunctor
     * methods are called or values corresponding to remote gps are stored.
     *
     * The input rGPDataMap map should not contain any key gps which are not included in creating the
     * GlobalPointerMapCommunicator used in this proxy (i.e. mrPointerCommunicator)
     *
     * This method can be called several times without severe additional computational cost
     * since this does not do any mpi communication. But in this case, remote gp values
     * will be overwritten, whereas local gp values will be applied using mrApplyFunctor.
     *
     * This does not have additional cost in serial run
     *
     * @see GlobalPointerMapCommunicator
     *
     * @param rGPDataMap        Input values map (key: GlobalPointer<TPointerDataType>, value: TValueDataType)
     */
    void Assign(const TGPDataMap& rGPDataMap)
    {
        // get gp vector for parallel omp run
        const auto& gps = TGPMapCommunicator::GetKeys(rGPDataMap);

        // running this in parallel assuming mrApplyFunctor is thread safe
        IndexPartition<int>(gps.size()).for_each([&](const int Index) {
            auto p_itr = rGPDataMap.find(gps[Index]);
            Assign(p_itr->first, p_itr->second);
        });
    }

    /**
     * @brief Assigns values with given list and values
     *
     * This method assigns local gps in rGPs using mrApplyFunctor instantly. Values
     * belonging to remote gps are stored in the mrNonLocalDataMap.
     *
     * In the case of serial job, this method calls AssignLocalData, where no checks are done
     * to ensure gps are local because all gps should be always local.
     *
     * In case of distributed job, this method calls AssignLocalAndRemoteData where checks are
     * performed to identify gps are local or non-local, based on that mrApplyFunctor
     * methods are called or values corresponding to remote gps are stored.
     *
     * The input rGPs map should not contain any key gps which are not included in creating the
     * GlobalPointerMapCommunicator used in this proxy (i.e. mrPointerCommunicator)
     *
     * This method can be called several times without severe additional computational cost
     * since this does not do any mpi communication. But in this case, remote gp values
     * will be overwritten, whereas local gp values will be applied using mrApplyFunctor.
     *
     * This does not have additional cost in serial run
     *
     * @see GlobalPointerMapCommunicator
     *
     * @param rGPs              List of GlobalPointers
     * @param rValues           List of values
     */
    void Assign(
        const TGPVector& rGPs,
        const std::vector<TValueDataType>& rValues)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rGPs.size() != rValues.size())
            << "Number of global pointers does not match with number of "
               "values. [ rGPs.size() = "
            << rGPs.size() << ", rValues.size() = " << rValues.size() << " ].\n";

        IndexPartition<int>(rGPs.size()).for_each([&](const int Index){
            Assign(rGPs(Index), rValues[Index]);
        });

        KRATOS_CATCH("");
    }

    /**
     * @brief Assigns value of the GlobalPointer
     *
     * This method assigns gp value by rValue using mrApplyFunctor in the case if
     * given rGlobalPointer is a local gp, otherwise mrNonLocalDataMap is assigned with
     * new value, which will be used in mpi communication.
     *
     * @see GlobalPointerMapCommunicator
     *
     * @param rGP               Global pointer of the destination
     * @param rValue            Value to be used in assigning
     */
    void Assign(
        const GlobalPointer<TPointerDataType>& rGlobalPointer,
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

    void (ProxyType::*mUpdateMethod)(const GlobalPointer<TPointerDataType>&, const TValueDataType&);
    void (ProxyType::*mNonLocalGPDataUpdateMethod)(TValueDataType&, const TValueDataType&);

    TGPNonLocalDataMap& mrNonLocalDataMap;
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
        const GlobalPointer<TPointerDataType>& rGlobalPointer,
        const TValueDataType& rValue)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rGlobalPointer.GetRank() != mCurrentRank)
            << "Using local global pointer update method with a non-local "
               "global pointer. [ MyPID = "
            << mCurrentRank
            << " GlobalPointerRank = " << rGlobalPointer.GetRank() << " ].\n";

        auto gp_pointer = rGlobalPointer;  // required to remove the const from rGlobalPointer
        mrApplyFunctor(*gp_pointer, rValue);

        KRATOS_CATCH("");
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
        const GlobalPointer<TPointerDataType>& rGlobalPointer,
        const TValueDataType& rValue)
    {
        KRATOS_TRY

        const int data_rank = rGlobalPointer.GetRank();

        if (data_rank == mCurrentRank) {
            AssignLocalData(rGlobalPointer, rValue);
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

            (this->*(this->mNonLocalGPDataUpdateMethod))(p_non_local_itr->second, rValue);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief This is used in assembly of values
     *
     * This method assembles non-local gp data in the non-local gp data map
     *
     * @param rNonLocalGPDataValue          Output value
     * @param rNewValue                     New input value
     */
    void NonLocalDataAssemblyMethod(
        TValueDataType& rNonLocalGPDataValue,
        const TValueDataType& rNewValue)
    {
        rNonLocalGPDataValue += rNewValue;
    }

    /**
     * @brief This is used in assigning of values
     *
     * This method assigns non-local gp data in the non-local gp data map
     *
     * @param rNonLocalGPDataValue          Output value
     * @param rNewValue                     New input value
     */
    void NonLocalDataAssignMethod(
        TValueDataType& rNonLocalGPDataValue,
        const TValueDataType& rNewValue)
    {
        rNonLocalGPDataValue = rNewValue;
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

    template<class TApplyFunctor>
    using ProxyType = ApplyProxy<TPointerDataType, TValueDataType, TApplyFunctor>;

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
     * This constructor will also compute communication scheduling
     * for given list of GPs.
     *
     * The rInitializationValue is used to initialize the gp map created
     * for GPs provided (for non local gps).
     *
     * @param rDataCommunicator         Data communicator
     * @param rGPVector                 Global pointers vector
     * @param rInitializationValue      Initialization value
     */
    GlobalPointerMapCommunicator(
        const DataCommunicator& rDataCommunicator,
        const TGPVector& rGPVector,
        const TValueDataType& rInitializationValue,
        const bool IsAssembly = true)
        : mrDataCommunicator(rDataCommunicator),
          mCurrentRank(rDataCommunicator.Rank()),
          mIsAssembly(IsAssembly)
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
     * TVectorFunctorType signature
     *      const DataCommunicator& -> GlobalPointersVector<TPointerDataType>
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
        TVectorFunctorType&& rVectorFunctor,
        const bool IsAssembly = true)
        : mrDataCommunicator(rDataCommunicator),
          mCurrentRank(rDataCommunicator.Rank()),
          mIsAssembly(IsAssembly)
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
        TVectorFunctorType&& rFunctor,
        const bool IsAssembly = true)
        : GlobalPointerMapCommunicator(
              ParallelEnvironment::GetDefaultDataCommunicator(),
              rInitializationValue,
              std::forward<TVectorFunctorType>(rFunctor),
              IsAssembly)
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
        if (IsDistributed()) {
            for (auto color : mColors) {
                if (color >= 0) {
                    const auto& received_gp_map = mrDataCommunicator.SendRecv(
                        mrNonLocalPointers[color], color, color);

                    // get gp vector for parallel omp run
                    const auto& gps = GetKeys(received_gp_map);

                    // running this in parallel assuming mrApplyFunctor is thread safe
                    IndexPartition<int>(gps.size()).for_each([&](const int Index) {
                        auto p_itr = received_gp_map.find(gps[Index]);
                        // it is safer to call the serial update method here
                        // because we should only get process's local gps when communication is done
                        rApplyProxy.AssignLocalData(p_itr->first, p_itr->second);
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
     * @tparam TApplyFunctor
     * @param rApplyFunctor
     * @return ApplyProxy<TPointerDataType, TValueDataType, TApplyFunctor>
     */
    template <class TApplyFunctor>
    ProxyType<TApplyFunctor> GetApplyProxy(TApplyFunctor&& rApplyFunctor)
    {
        return ProxyType<TApplyFunctor>(
            std::forward<TApplyFunctor>(rApplyFunctor), mrNonLocalPointers, *this);
    }

    /**
     * @brief Returns jobs parallel status
     *
     */
    bool IsDistributed() const
    {
        return mrDataCommunicator.IsDistributed();
    }

    int GetMyPID() const
    {
        return mCurrentRank;
    }

    bool IsAssembly() const
    {
        return mIsAssembly;
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
        // this method is called always at the constructor.
        // this should only be called in MPI
        for (const auto& r_gp : rGPVector.GetContainer()) {
            if (r_gp.GetRank() != mCurrentRank) {
                    auto& rank_local_map = mrNonLocalPointers[r_gp.GetRank()];
                    rank_local_map[r_gp] = rInitializationValue;
            }
        }
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
    const bool mIsAssembly;

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


