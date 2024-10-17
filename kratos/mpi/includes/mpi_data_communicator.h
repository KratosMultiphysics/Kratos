//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <mpi.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SYNC_SHAPE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SYNC_SHAPE_INTERFACE_FOR_TYPE(...)         \
bool SynchronizeShape(__VA_ARGS__&) const override;                                     \
bool SynchronizeShape(                                                                  \
    const __VA_ARGS__& rSendValue, const int SendDestination, const int SendTag,        \
    __VA_ARGS__& rRecvValue, const int RecvSource, const int RecvTag) const override;   \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(...)                                                         \
__VA_ARGS__ Sum(const __VA_ARGS__& rLocalValue, const int Root) const override;                                                     \
std::vector<__VA_ARGS__> Sum(const std::vector<__VA_ARGS__>& rLocalValues, const int Root) const override;                          \
void Sum(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues, const int Root) const override;     \
__VA_ARGS__ Min(const __VA_ARGS__& rLocalValue, const int Root) const override;                                                     \
std::vector<__VA_ARGS__> Min(const std::vector<__VA_ARGS__>& rLocalValues, const int Root) const override;                          \
void Min(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues, const int Root) const override;     \
__VA_ARGS__ Max(const __VA_ARGS__& rLocalValue, const int Root) const override;                                                     \
std::vector<__VA_ARGS__> Max(const std::vector<__VA_ARGS__>& rLocalValues, const int Root) const override;                          \
void Max(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues, const int Root) const override;     \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(...)                                      \
__VA_ARGS__ SumAll(const __VA_ARGS__& rLocalValue) const override;                                                  \
std::vector<__VA_ARGS__> SumAll(const std::vector<__VA_ARGS__>& rLocalValues) const override;                       \
void SumAll(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const override;  \
__VA_ARGS__ MinAll(const __VA_ARGS__& rLocalValue) const override;                                                  \
std::vector<__VA_ARGS__> MinAll(const std::vector<__VA_ARGS__>& rLocalValues) const override;                       \
void MinAll(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const override;  \
__VA_ARGS__ MaxAll(const __VA_ARGS__& rLocalValue) const override;                                                  \
std::vector<__VA_ARGS__> MaxAll(const std::vector<__VA_ARGS__>& rLocalValues) const override;                       \
void MaxAll(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const override;  \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE(...)                                  \
std::pair<__VA_ARGS__, int> MinLocAll(const __VA_ARGS__& rLocalValue) const override;                               \
std::pair<__VA_ARGS__, int> MaxLocAll(const __VA_ARGS__& rLocalValue) const override;
#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(...)                                            \
__VA_ARGS__ ScanSum(const __VA_ARGS__& rLocalValue) const override;                                                     \
std::vector<__VA_ARGS__> ScanSum(const std::vector<__VA_ARGS__>& rLocalValues) const override;                          \
void ScanSum(const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const override;     \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE(...)                       \
__VA_ARGS__ SendRecvImpl(                                                                           \
    const __VA_ARGS__& SendValue, const int SendDestination, const int SendTag,                     \
    const int RecvSource, const int RecvTag) const override;                                        \
std::vector<__VA_ARGS__> SendRecvImpl(const std::vector<__VA_ARGS__>& rSendValues,                  \
    const int SendDestination, const int SendTag,                                                   \
    const int RecvSource, const int RecvTag) const override;                                        \
void SendRecvImpl(                                                                                  \
    const __VA_ARGS__& SendValue, const int SendDestination, const int SendTag,                     \
    __VA_ARGS__& RecvValue, const int RecvSource, const int RecvTag) const override;                \
void SendRecvImpl(                                                                                  \
    const std::vector<__VA_ARGS__>& rSendValues, const int SendDestination, const int SendTag,      \
    std::vector<__VA_ARGS__>& rRecvValues, const int RecvSource, const int RecvTag) const override; \
void SendImpl(const __VA_ARGS__& rSendValues,                                                       \
    const int SendDestination, const int SendTag = 0) const override;                               \
void SendImpl(const std::vector<__VA_ARGS__>& rSendValues,                                          \
    const int SendDestination, const int SendTag = 0) const override;                               \
void RecvImpl(__VA_ARGS__& rRecvValues,                                                             \
    const int RecvSource, const int RecvTag = 0) const override;                                    \
void RecvImpl(std::vector<__VA_ARGS__>& rRecvValues,                                                \
    const int RecvSource, const int RecvTag = 0) const override;                                    \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE(...)                  \
void BroadcastImpl(__VA_ARGS__& rBuffer, const int SourceRank) const override;                  \
void BroadcastImpl(std::vector<__VA_ARGS__>& rBuffer, const int SourceRank) const override;     \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE(...)                        \
std::vector<__VA_ARGS__> Scatter(                                                                   \
    const std::vector<__VA_ARGS__>& rSendValues, const int SourceRank) const override;              \
void Scatter(                                                                                       \
    const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues,             \
    const int SourceRank) const override;                                                           \
std::vector<__VA_ARGS__> Scatterv(                                                                  \
    const std::vector<std::vector<__VA_ARGS__>>& rSendValues, const int SourceRank) const override; \
void Scatterv(                                                                                      \
    const std::vector<__VA_ARGS__>& rSendValues,                                                    \
    const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,                      \
    std::vector<__VA_ARGS__>& rRecvValues, const int SourceRank) const override;                    \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE(...)                                             \
std::vector<__VA_ARGS__> Gather(const std::vector<__VA_ARGS__>& rSendValues, const int DestinationRank) const override; \
void Gather(                                                                                                            \
    const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues,                                 \
    const int DestinationRank) const override;                                                                          \
std::vector<std::vector<__VA_ARGS__>> Gatherv(                                                                          \
    const std::vector<__VA_ARGS__>& rSendValues, const int DestinationRank) const override;                             \
void Gatherv(const std::vector<__VA_ARGS__>& rSendValues,                                                               \
        std::vector<__VA_ARGS__>& rRecvValues,                                                                          \
        const std::vector<int>& rRecvCounts,                                                                            \
        const std::vector<int>& rRecvOffsets,                                                                           \
        const int DestinationRank) const override;                                                                      \
std::vector<__VA_ARGS__> AllGather(const std::vector<__VA_ARGS__>& rSendValues) const override;                         \
void AllGather(const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues) const override;      \
std::vector<std::vector<__VA_ARGS__>> AllGatherv(const std::vector<__VA_ARGS__>& rSendValues) const  override;          \
void AllGatherv(const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues,                     \
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const override;
#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(...)     \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(__VA_ARGS__)     \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(__VA_ARGS__)  \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(__VA_ARGS__)    \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE(__VA_ARGS__)    \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE(__VA_ARGS__)     \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SYNC_SHAPE_INTERFACE_FOR_TYPE(__VA_ARGS__) \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(...)   \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE(__VA_ARGS__)  \
KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE(__VA_ARGS__) \

#endif

namespace Kratos
{
///@addtogroup Kratos MPI Core
///@{

///@name Kratos Classes
///@{

/// Wrapper for common MPI calls within Kratos.
/** This class is designed to isolate the Kratos core and applications from direct calls to MPI routines.
 *
 *  For function operating on std::vectors, no effort is made to resize inconsistent vectors.
 *  The sizes of the sending and receiving buffers will only be checked if the code is compiled in Debug
 *  mode. In that case, a meaningful error explaining the inconsistency will be produced.
 *  This is done for efficiency (size checks can force multi-stage communications).
 *
 *  @see DataCommunicator in the KratosCore for the full interface and a serial do-nothing implementation.
 */
class KRATOS_API(KRATOS_MPI_CORE) MPIDataCommunicator: public DataCommunicator
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPIDataCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(MPIDataCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor accepting an MPI_Comm object.
    explicit MPIDataCommunicator(MPI_Comm MPIComm);

    /// Destructor.
    ~MPIDataCommunicator() override;

    ///@}
    ///@name Operations
    ///@{

    /// Create a new MPIDataCommunicator using the provided MPI_Comm object.
    /** The new MPIDataCommunicator instance is returned as a unique pointer,
     *  since it is responsible for managing the lifetime of the underlying MPI_Comm,
     *  and in particular calling MPI_Comm_free once it goes out of scope
     *  (this is only required/done if Comm is not one of the predefined MPI_COMM types).
     */
    static MPIDataCommunicator::UniquePointer Create(MPI_Comm MPIComm);

    void Barrier() const override;

    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(char)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(unsigned int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(long unsigned int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(double)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 3>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 4>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 6>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 9>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(Vector)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(Matrix)

    // MinLoc and MaxLoc AllReduce operations
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE(char)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE(int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE(unsigned int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE(long unsigned int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_LOC_INTERFACE_FOR_TYPE(double)

    // Reduce operations

    bool AndReduce(
        const bool Value,
        const int Root) const override;

    Kratos::Flags AndReduce(
        const Kratos::Flags Values,
        const Kratos::Flags Mask,
        const int Root) const override;

    bool OrReduce(
        const bool Value,
        const int Root) const override;

    Kratos::Flags OrReduce(
        const Kratos::Flags Values,
        const Kratos::Flags Mask,
        const int Root) const override;

    // Allreduce operations

    bool AndReduceAll(const bool Value) const override;

    Kratos::Flags AndReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const override;

    bool OrReduceAll(const bool Value) const override;

    Kratos::Flags OrReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const override;

    ///@}
    ///@name Access
    ///@{

    /// Get the underlying MPI_Comm instance
    /** @note This method does not exist in the base class
     *  as it would introduce a dependency to MPI in the Kratos core.
     */
    static MPI_Comm GetMPICommunicator(const DataCommunicator& rDataCommunicator);

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the parallel rank for this DataCommunicator.
     * @details This function serves as a wrapper for MPI_Comm_rank.
     * @return The parallel rank of the current process.
     */
    int Rank() const override;

    /**
     * @brief Get the parallel size of this DataCommunicator.
     * @details This function serves as a wrapper for MPI_Comm_size.
     * @return The parallel size of the communicator.
     */
    int Size() const override;

    /**
     * @brief Check whether this DataCommunicator is aware of parallelism.
     * @return True if the DataCommunicator is distributed, otherwise false.
     */
    bool IsDistributed() const override;

    /**
     * @brief Check whether this DataCommunicator involves the current rank.
     * @details In MPI, if the rank is not involved in communication, the communicator is MPI_COMM_NULL and is not a valid argument for most MPI calls.
     * @return True if the DataCommunicator is defined on the current rank, otherwise false.
     */
    bool IsDefinedOnThisRank() const override;

    /**
     * @brief Check whether this DataCommunicator is MPI_COMM_NULL for the current rank.
     * @details In MPI, if the rank is not involved in communication, the communicator is MPI_COMM_NULL and is not a valid argument for most MPI calls.
     * @return True if the DataCommunicator is MPI_COMM_NULL, otherwise false.
     */
    bool IsNullOnThisRank() const override;

    /**
     * @brief Get a sub-data communicator.
     * @details This function returns a sub-data communicator based on the provided ranks and a new communicator name.
     * @param rRanks               The ranks to include in the sub-communicator.
     * @param rNewCommunicatorName The name of the new sub-communicator.
     * @return The sub-data communicator.
     */
    const DataCommunicator& GetSubDataCommunicator(
        const std::vector<int>& rRanks,
        const std::string& rNewCommunicatorName
        ) const override;

    ///@}
    ///@name Helper functions for error checking in MPI
    ///@{

    bool BroadcastErrorIfTrue(bool Condition, const int SourceRank) const override;

    bool BroadcastErrorIfFalse(bool Condition, const int SourceRank) const override;

    bool ErrorIfTrueOnAnyRank(bool Condition) const override;

    bool ErrorIfFalseOnAnyRank(bool Condition) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override;

    ///@}

  protected:

    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(char)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(unsigned int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(long unsigned int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(double)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 3>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 4>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 6>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 9>)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(Vector)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(Matrix)

    // Broadcast operations

    void BroadcastImpl(std::string& rBroadcastValues, const int SourceRank) const override;

    // Sendrecv operations

    std::string SendRecvImpl(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        const int RecvSource, const int RecvTag) const override;

    void SendRecvImpl(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        std::string& rRecvValues, const int RecvSource, const int RecvTag) const override;

    void SendImpl(const std::string& rSendValues, const int SendDestination, const int SendTag = 0) const override;

    void RecvImpl(std::string& rRecvValues, const int RecvSource, const int RecvTag = 0) const override;

  private:
    ///@name Member Variables
    ///@{

    MPI_Comm mComm;

    ///@}
    ///@name Operations
    ///@{

    void CheckMPIErrorCode(const int ierr, const std::string& MPICallName) const;

    template<class TDataType> bool SynchronizeShapeDetail(TDataType& rValue) const;

    template<class TDataType> bool SynchronizeShapeDetail(
        const TDataType& rSendValue,
        const int SendDestination,
        const int SendTag,
        TDataType& rRecvValue,
        const int RecvSource,
        const int RecvTag) const;

    template<class TDataType> void ReduceDetail(
        const TDataType& rLocalValues,
        TDataType& rReducedValues,
        MPI_Op Operation,
        const int Root) const;

    template<class TDataType> TDataType ReduceDetail(
        const TDataType& rLocalValues,
        MPI_Op Operation,
        const int Root) const;

    template<class TDataType> std::vector<TDataType> ReduceDetailVector(
        const std::vector<TDataType>& rLocalValues,
        MPI_Op Operation,
        const int Root) const;

    template<class TDataType> void AllReduceDetail(
        const TDataType& rLocalValues,
        TDataType& rReducedValues,
        MPI_Op Operation) const;

    template<class TDataType> TDataType AllReduceDetail(
        const TDataType& rLocalValues, MPI_Op Operation) const;

    template<class TDataType> std::vector<TDataType> AllReduceDetailVector(
        const std::vector<TDataType>& rLocalValues,
        MPI_Op Operation) const;

    /**
    * @brief Performs an AllReduce operation with location information (the partition where the reduced value was found).
    * @details This function performs an AllReduce operation on a pair of data and an integer location using the specified MPI operation. The AllReduce operation combines the data from all processes and stores the result in the pair's first element. The location information (integer) is the partition where the reduced value was found.
    * @tparam TDataType The data type of the pair's first element.
    * @param rLocalValues A pair containing the local data and location information to be reduced.
    * @param Operation The MPI operation to use for the reduction.
    * @return A pair where the first element contains the result of the AllReduce operation, and the second element is the partition where the reduced value was found.
    */
    template<class TDataType>
    std::pair<TDataType, int> AllReduceDetailWithLocation(
        const std::pair<TDataType, int>& rLocalValues,
        MPI_Op Operation
        ) const;

    template<class TDataType> void ScanDetail(
        const TDataType& rLocalValues,
        TDataType& rReducedValues,
        MPI_Op Operation) const;

    template<class TDataType> TDataType ScanDetail(
        const TDataType rLocalValues,
        MPI_Op Operation) const;

    template<class TDataType> std::vector<TDataType> ScanDetail(
        const std::vector<TDataType>& rLocalValues,
        MPI_Op Operation) const;

    template<class TDataType> void SendRecvDetail(
        const TDataType& rSendMessage, const int SendDestination, const int SendTag,
        TDataType& rRecvMessage, const int RecvSource, const int RecvTag) const;

    template<class TDataType> TDataType SendRecvDetail(
        const TDataType& rSendMessage,
        const int SendDestination, const int SendTag,
        const int RecvSource, const int RecvTag) const;

    template<class TDataType> std::vector<TDataType> SendRecvDetail(
        const std::vector<TDataType>& rSendMessage,
        const int SendDestination, const int SendTag,
        const int RecvSource, const int RecvTag) const;

    template<class TDataType> void SendDetail(
        const TDataType& rSendValues, const int SendDestination, const int SendTag) const;

    template<class TDataType> void RecvDetail(
        TDataType& rRecvValues, const int RecvSource, const int RecvTag) const;

    template<class TDataType> void BroadcastDetail(
        TDataType& rBuffer, const int SourceRank) const;

    template<class TSendDataType, class TRecvDataType> void ScatterDetail(
        const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int SourceRank) const;

    template<class TDataType> std::vector<TDataType> ScatterDetail(
        const std::vector<TDataType>& rSendValues, const int SourceRank) const;

    template<class TDataType> void ScattervDetail(
        const TDataType& rSendValues,
        const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const;

    template<class TDataType> std::vector<TDataType> ScattervDetail(
        const std::vector<std::vector<TDataType>>& rSendValues,const int SourceRank) const;

    template<class TSendDataType, class TRecvDataType> void GatherDetail(
        const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int RecvRank) const;

    template<class TDataType> std::vector<TDataType> GatherDetail(
        const std::vector<TDataType>& rSendValues, const int DestinationRank) const;

    template<class TDataType> void GathervDetail(
        const TDataType& rSendValues, TDataType& rRecvValues,
        const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
        const int RecvRank) const;

    template<class TDataType> std::vector<std::vector<TDataType>>
    GathervDetail(const std::vector<TDataType>& rSendValues, const int DestinationRank) const;

    template<class TDataType> void AllGatherDetail(
        const TDataType& rSendValues, TDataType& rRecvValues) const;

    template<class TDataType> std::vector<TDataType> AllGatherDetail(
        const std::vector<TDataType>& rSendValues) const;

    template<class TDataType>
    void AllGathervDetail(
        const TDataType& rSendValues, TDataType& rRecvValues,
        const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const;

    template<class TDataType>
    std::vector<std::vector<TDataType>> AllGathervDetail(
        const std::vector<TDataType>& rSendValues) const;

    bool IsEqualOnAllRanks(const int LocalValue) const;

    bool IsValidRank(const int Rank) const;

    template<class TDataType> void ValidateScattervInput(
        const TDataType& rSendValues,
        const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const;

    template<class TDataType> void ValidateGathervInput(
        const TDataType& rSendValues, TDataType& rRecvValues,
        const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
        const int RecvRank) const;

    template<class TDataType>
    void ValidateAllGathervInput(
        const TDataType& rSendValues, TDataType& rRecvValues,
        const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const;

    template<class TDataType> void PrepareScattervBuffers(
        const std::vector<std::vector<TDataType>>& rInputMessage,
        std::vector<TDataType>& rScattervMessage,
        std::vector<int>& rMessageLengths,
        std::vector<int>& rMessageDistances,
        std::vector<TDataType>& rResult,
        const int SourceRank) const;

    template<class TDataType> void PrepareGathervBuffers(
        const std::vector<TDataType>& rGathervInput,
        std::vector<TDataType>& rGathervMessage,
        std::vector<int>& rMessageLengths,
        std::vector<int>& rMessageDistances,
        const int DestinationRank) const;

    template<class TDataType>
    void PrepareAllGathervBuffers(
        const std::vector<TDataType>& rGathervInput,
        std::vector<TDataType>& rGathervMessage,
        std::vector<int>& rMessageLengths,
        std::vector<int>& rMessageDistances) const;

    template<class TDataType> void PrepareGathervReturn(
        const std::vector<TDataType>& rGathervMessage,
        const std::vector<int>& rMessageLengths,
        const std::vector<int>& rMessageDistances,
        std::vector<std::vector<TDataType>>& rOutputMessage,
        const int DestinationRank) const;

    template<class TDataType>
    void PrepareAllGathervReturn(
        const std::vector<TDataType>& rGathervMessage,
        const std::vector<int>& rMessageLengths,
        const std::vector<int>& rMessageDistances,
        std::vector<std::vector<TDataType>>& rOutputMessage) const;

    template<class TValue> inline MPI_Datatype MPIDatatype(const TValue&) const;

    template<class TContainer> inline void* MPIBuffer(TContainer& rValues) const;

    template<class TContainer> inline const void* MPIBuffer(const TContainer& rValues) const;

    template<class TContainer> inline int MPIMessageSize(const TContainer& rValues) const;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Copy constructor.
    MPIDataCommunicator(MPIDataCommunicator const &rOther) = delete;

    /// Assignment operator.
    MPIDataCommunicator &operator=(MPIDataCommunicator const &rOther) = delete;

    ///@}

}; // Class MPIDataCommunicator

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                MPIDataCommunicator &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const MPIDataCommunicator &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE
