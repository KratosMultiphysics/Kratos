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

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED
#define KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <mpi.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(type)                                      \
type Sum(const type rLocalValue, const int Root) const override;                                                  \
std::vector<type> Sum(const std::vector<type>& rLocalValues, const int Root) const override;                      \
void Sum(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const override; \
type Min(const type rLocalValue, const int Root) const override;                                                  \
std::vector<type> Min(const std::vector<type>& rLocalValues, const int Root) const override;                      \
void Min(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const override; \
type Max(const type rLocalValue, const int Root) const override;                                                  \
std::vector<type> Max(const std::vector<type>& rLocalValues, const int Root) const override;                      \
void Max(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const override; \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(type)                      \
type SumAll(const type rLocalValue) const override;                                                  \
std::vector<type> SumAll(const std::vector<type>& rLocalValues) const override;                      \
void SumAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const override; \
type MinAll(const type rLocalValue) const override;                                                  \
std::vector<type> MinAll(const std::vector<type>& rLocalValues) const override;                      \
void MinAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const override; \
type MaxAll(const type rLocalValue) const override;                                                  \
std::vector<type> MaxAll(const std::vector<type>& rLocalValues) const override;                      \
void MaxAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const override; \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(type)                         \
type ScanSum(const type rLocalValue) const override;                                                  \
std::vector<type> ScanSum(const std::vector<type>& rLocalValues) const override;                      \
void ScanSum(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const override; \

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
class MPIDataCommunicator: public DataCommunicator
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

    DataCommunicator::UniquePointer Clone() const override;

    void Barrier() const override;

    // Reduce operations

    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(double)

    array_1d<double,3> Sum(const array_1d<double,3>& rLocalValue, const int Root) const override;

    array_1d<double,3> Min(const array_1d<double,3>& rLocalValue, const int Root) const override;

    array_1d<double,3> Max(const array_1d<double,3>& rLocalValue, const int Root) const override;

    Kratos::Flags AndReduce(
        const Kratos::Flags Values,
        const Kratos::Flags Mask,
        const int Root) const override;

    Kratos::Flags OrReduce(
        const Kratos::Flags Values,
        const Kratos::Flags Mask,
        const int Root) const override;

    // Allreduce operations

    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(double)

    array_1d<double,3> SumAll(const array_1d<double,3>& rLocalValue) const override;

    array_1d<double,3> MinAll(const array_1d<double,3>& rLocalValue) const override;

    array_1d<double,3> MaxAll(const array_1d<double,3>& rLocalValue) const override;

    Kratos::Flags AndReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const override;

    Kratos::Flags OrReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const override;

    // Scan operations

    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(int)
    KRATOS_MPI_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(double)

    // Sendrecv operations

    std::vector<int> SendRecv(
        const std::vector<int>& rSendValues,
        const int SendDestination,
        const int RecvSource) const override;

    std::vector<double> SendRecv(
        const std::vector<double>& rSendValues,
        const int SendDestination,
        const int RecvSource) const override;

    std::string SendRecv(
        const std::string& rSendValues,
        const int SendDestination,
        const int RecvSource) const override;

    void SendRecv(
        const std::vector<int>& rSendValues, const int SendDestination, const int SendTag,
        std::vector<int>& rRecvValues, const int RecvSource, const int RecvTag) const override;

    void SendRecv(
        const std::vector<double>& rSendValues, const int SendDestination, const int SendTag,
        std::vector<double>& rRecvValues, const int RecvSource, const int RecvTag) const override;

    void SendRecv(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        std::string& rRecvValues, const int RecvSource, const int RecvTag) const override;

    // Broadcast

    void Broadcast(
        int& rBuffer,
        const int SourceRank) const override;

    void Broadcast(
        double& rBuffer,
        const int SourceRank) const override;

    void Broadcast(
        std::vector<int>& rBuffer,
        const int SourceRank) const override;

    void Broadcast(
        std::vector<double>& rBuffer,
        const int SourceRank) const override;

    // Scatter operations

    std::vector<int> Scatter(
        const std::vector<int>& rSendValues,
        const int SourceRank) const override;

    std::vector<double> Scatter(
        const std::vector<double>& rSendValues,
        const int SourceRank) const override;

    void Scatter(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const int SourceRank) const override;

    void Scatter(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const int SourceRank) const override;

    // Scatterv operations

    std::vector<int> Scatterv(
        const std::vector<std::vector<int>>& rSendValues,
        const int SourceRank) const override;

    std::vector<double> Scatterv(
        const std::vector<std::vector<double>>& rSendValues,
        const int SourceRank) const override;

    void Scatterv(
        const std::vector<int>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<int>& rRecvValues,
        const int SourceRank) const override;

    void Scatterv(
        const std::vector<double>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<double>& rRecvValues,
        const int SourceRank) const override;

    // Gather operations

    std::vector<int> Gather(
        const std::vector<int>& rSendValues,
        const int DestinationRank) const override;

    std::vector<double> Gather(
        const std::vector<double>& rSendValues,
        const int DestinationRank) const override;

    void Gather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const int DestinationRank) const override;

    void Gather(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const int DestinationRank) const override;

    // Gatherv operations

    std::vector<std::vector<int>> Gatherv(
        const std::vector<int>& rSendValues,
        const int DestinationRank) const override;

    std::vector<std::vector<double>> Gatherv(
        const std::vector<double>& rSendValues,
        const int DestinationRank) const override;

    void Gatherv(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const int DestinationRank) const override;

    void Gatherv(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const int DestinationRank) const override;

    // Allgather operations

    std::vector<int> AllGather(
        const std::vector<int>& rSendValues) const override;

    std::vector<double> AllGather(
        const std::vector<double>& rSendValues) const override;

    void AllGather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues) const override;

    void AllGather(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues) const override;

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

    int Rank() const override;

    int Size() const override;

    bool IsDistributed() const override;

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

  private:
    ///@name Member Variables
    ///@{

    MPI_Comm mComm;

    ///@}
    ///@name Operations
    ///@{

    void CheckMPIErrorCode(const int ierr, const std::string& MPICallName) const;

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

    template<class TDataType> void BroadcastDetail(
        TDataType& rBuffer, const int SourceRank) const;

    template<class TSendDataType, class TRecvDataType> void ScatterDetail(
        const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int SourceRank) const;

    template<class TDataType> void ScattervDetail(
        const TDataType& rSendValues,
        const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const;

    template<class TSendDataType, class TRecvDataType> void GatherDetail(
        const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int RecvRank) const;

    template<class TDataType> void GathervDetail(
        const TDataType& rSendValues, TDataType& rRecvValues,
        const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
        const int RecvRank) const;

    template<class TDataType> void AllGatherDetail(
        const TDataType& rSendValues, TDataType& rRecvValues) const;

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

    template<class TDataType> void PrepareGathervReturn(
        const std::vector<TDataType>& rGathervMessage,
        const std::vector<int>& rMessageLengths,
        const std::vector<int>& rMessageDistances,
        std::vector<std::vector<TDataType>>& rOutputMessage,
        const int DestinationRank) const;

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

#endif // KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED  defined