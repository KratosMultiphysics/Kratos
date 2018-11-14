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

namespace Internals {

template<class TValue> inline MPI_Datatype MPIDatatype(const TValue&);

}

///@name Kratos Classes
///@{

/// Wrapper for common MPI calls within Kratos.
/** This class is designed to isolate the Kratos core and applications from direct calls to MPI routines.
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

    void SumAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const override;

    void SumAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const override;

    int MinAll(const int rLocalValue) const override;

    double MinAll(const double rLocalValue) const override;

    array_1d<double,3> MinAll(const array_1d<double,3>& rLocalValue) const override;

    void MinAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const override;

    void MinAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const override;

    int MaxAll(const int rLocalValue) const override;

    double MaxAll(const double rLocalValue) const override;

    array_1d<double,3> MaxAll(const array_1d<double,3>& rLocalValue) const override;

    void MaxAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const override;

    void MaxAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const override;

    // Scan operations

    int ScanSum(const int rLocalValue) const override;

    double ScanSum(const double rLocalValue) const override;

    void ScanSum(const std::vector<int>& rLocalValues, std::vector<int>& rPartialSums) const override;

    void ScanSum(const std::vector<double>& rLocalValues, std::vector<double>& rPartialSums) const override;

    // Sendrecv operations

    void SendRecv(
        const std::vector<int>& rSendValues, const unsigned int SendDestination,
        std::vector<int>& rRecvValues, const unsigned int RecvSource) const override;

    void SendRecv(
        const std::vector<double>& rSendValues, const unsigned int SendDestination,
        std::vector<double>& rRecvValues, const unsigned int RecvSource) const override;

    // Broadcast

    void Broadcast(
        int& rBuffer,
        const unsigned int SourceRank) const override;

    void Broadcast(
        double& rBuffer,
        const unsigned int SourceRank) const override;

    void Broadcast(
        std::vector<int>& rBuffer,
        const unsigned int SourceRank) const override;

    void Broadcast(
        std::vector<double>& rBuffer,
        const unsigned int SourceRank) const override;

    // Scatter operations

    void Scatter(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const unsigned int SourceRank) const override;

    void Scatter(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const unsigned int SourceRank) const override;

    void Scatterv(
        const std::vector<int>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<int>& rRecvValues,
        const unsigned int SourceRank) const override;

    void Scatterv(
        const std::vector<double>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<double>& rRecvValues,
        const unsigned int SourceRank) const override;

    // Gather operations

    void Gather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const unsigned int DestinationRank) const override;

    void Gather(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const unsigned int DestinationRank) const override;

    void Gatherv(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const unsigned int DestinationRank) const override;

    void Gatherv(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const unsigned int DestinationRank) const override;

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
