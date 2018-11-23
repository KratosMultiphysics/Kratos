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
    MPIDataCommunicator(MPI_Comm MPIComm);

    /// Destructor.
    ~MPIDataCommunicator() override;

    ///@}
    ///@name Operations
    ///@{

    DataCommunicator::UniquePointer Clone() const override;

    void Barrier() const override;

    // Reduce operations

    int Sum(const int rLocalValue, const int Root) const override;

    double Sum(const double rLocalValue, const int Root) const override;

    array_1d<double,3> Sum(const array_1d<double,3>& rLocalValue, const int Root) const override;

    std::vector<int> Sum(const std::vector<int>& rLocalValues, const int Root) const override;

    std::vector<double> Sum(const std::vector<double>& rLocalValues, const int Root) const override;

    void Sum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const override;

    void Sum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const override;

    int Min(const int rLocalValue, const int Root) const override;

    double Min(const double rLocalValue, const int Root) const override;

    array_1d<double,3> Min(const array_1d<double,3>& rLocalValue, const int Root) const override;

    std::vector<int> Min(const std::vector<int>& rLocalValues, const int Root) const override;

    std::vector<double> Min(const std::vector<double>& rLocalValues, const int Root) const override;

    void Min(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const override;

    void Min(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const override;

    int Max(const int rLocalValue, const int Root) const override;

    double Max(const double rLocalValue, const int Root) const override;

    array_1d<double,3> Max(const array_1d<double,3>& rLocalValue, const int Root) const override;

    std::vector<int> Max(const std::vector<int>& rLocalValues, const int Root) const override;

    std::vector<double> Max(const std::vector<double>& rLocalValues, const int Root) const override;

    void Max(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const override;

    void Max(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const override;

    // Allreduce operations

    int SumAll(const int rLocalValue) const override;

    double SumAll(const double rLocalValue) const override;

    array_1d<double,3> SumAll(const array_1d<double,3>& rLocalValue) const override;

    std::vector<int> SumAll(const std::vector<int>& rLocalValue) const override;

    std::vector<double> SumAll(const std::vector<double>& rLocalValue) const override;

    void SumAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const override;

    void SumAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const override;

    int MinAll(const int rLocalValue) const override;

    double MinAll(const double rLocalValue) const override;

    array_1d<double,3> MinAll(const array_1d<double,3>& rLocalValue) const override;

    std::vector<int> MinAll(const std::vector<int>& rLocalValues) const override;

    std::vector<double> MinAll(const std::vector<double>& rLocalValues) const override;

    void MinAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const override;

    void MinAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const override;

    int MaxAll(const int rLocalValue) const override;

    double MaxAll(const double rLocalValue) const override;

    array_1d<double,3> MaxAll(const array_1d<double,3>& rLocalValue) const override;

    std::vector<int> MaxAll(const std::vector<int>& rLocalValues) const override;

    std::vector<double> MaxAll(const std::vector<double>& rLocalValues) const override;

    void MaxAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const override;

    void MaxAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const override;

    // Scan operations

    int ScanSum(const int rLocalValue) const override;

    double ScanSum(const double rLocalValue) const override;

    std::vector<int> ScanSum(const std::vector<int>& rLocalValues) const override;

    std::vector<double> ScanSum(const std::vector<double>& rLocalValues) const override;

    void ScanSum(const std::vector<int>& rLocalValues, std::vector<int>& rPartialSums) const override;

    void ScanSum(const std::vector<double>& rLocalValues, std::vector<double>& rPartialSums) const override;

    // Sendrecv operations

    std::vector<int> SendRecv(
        const std::vector<int>& rSendValues,
        const int SendDestination,
        const int RecvSource) const override;

    std::vector<double> SendRecv(
        const std::vector<double>& rSendValues,
        const int SendDestination,
        const int RecvSource) const override;

    void SendRecv(
        const std::vector<int>& rSendValues, const int SendDestination,
        std::vector<int>& rRecvValues, const int RecvSource) const override;

    void SendRecv(
        const std::vector<double>& rSendValues, const int SendDestination,
        std::vector<double>& rRecvValues, const int RecvSource) const override;

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
    MPI_Comm GetMPICommunicator() const;

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

    void CheckMPIErrorCode(const int ierr, const std::string MPICallName) const;

    template<class TDataType> void ReduceDetail(
        const TDataType& rLocalValues,
        TDataType& rReducedValues,
        MPI_Op Operation,
        const int Root) const;

    template<class TDataType> void AllReduceDetail(
        const TDataType& rLocalValues,
        TDataType& rReducedValues,
        MPI_Op Operation) const;

    template<class TDataType> void ScanDetail(
        const TDataType& rLocalValues,
        TDataType& rReducedValues,
        MPI_Op Operation) const;

    template<class TDataType> void SendRecvDetail(
        const TDataType& rSendMessage, const int SendDestination,
        TDataType& rRecvMessage, const int RecvSource) const;

    template<class TDataType> void BroadcastDetail(
        TDataType& rBuffer, const int SourceRank) const;

    template<class TDataType> void ScatterDetail(
        const TDataType& rSendValues, TDataType& rRecvValues, const int SourceRank) const;

    template<class TDataType> void ScattervDetail(
        const TDataType& rSendValues,
        const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const;

    template<class TDataType> void GatherDetail(
        const TDataType& rSendValues, TDataType& rRecvValues, const int RecvRank) const;

    template<class TDataType> void GathervDetail(
        const TDataType& rSendValues, TDataType& rRecvValues,
        const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
        const int RecvRank) const;

    template<class TDataType> void AllGatherDetail(
        const TDataType& rSendValues, TDataType& rRecvValues) const;

    bool IsEqualOnAllRanks(const int LocalValue) const;

    bool IsValidRank(const int Rank) const;

    template<class TDataType> void ValidateSendRecvInput(
        const TDataType& rSendMessage, const int SendDestination,
        TDataType& rRecvMessage, const int RecvSource) const;

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
        const int SourceRank) const;

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

#endif // KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED  defined
