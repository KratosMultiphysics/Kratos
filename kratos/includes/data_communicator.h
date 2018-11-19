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

#ifndef KRATOS_DATA_COMMUNICATOR_H_INCLUDED
#define KRATOS_DATA_COMMUNICATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Serial (do-nothing) version of a wrapper class for MPI communication.
/** @see MPIDataCommunicator for a working distributed memory implementation.
  */
class DataCommunicator
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(DataCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataCommunicator() {}

    /// Destructor.
    virtual ~DataCommunicator() {}

    ///@}
    ///@name Operations
    ///@{

    virtual DataCommunicator::UniquePointer Clone() const
    {
        return Kratos::make_unique<DataCommunicator>();
    }

    virtual void Barrier() const {}

    // Reduce operations

    virtual int Sum(const int rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual double Sum(const double rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> Sum(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual void Sum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
    {}

    virtual void Sum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
    {}

    virtual int Min(const int rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual double Min(const double rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> Min(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual void Min(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
    {}

    virtual void Min(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
    {}

    virtual int Max(const int rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual double Max(const double rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> Max(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual void Max(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
    {}

    virtual void Max(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
    {}

    // Allreduce operations

    virtual int SumAll(const int rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double SumAll(const double rLocalValue) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> SumAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual void SumAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
    {}

    virtual void SumAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
    {}

    virtual int MinAll(const int rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double MinAll(const double rLocalValue) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> MinAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual void MinAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
    {}

    virtual void MinAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
    {}

    virtual int MaxAll(const int rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double MaxAll(const double rLocalValue) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> MaxAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual void MaxAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
    {}

    virtual void MaxAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
    {}

    // Scan operations

    virtual int ScanSum(const int rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double ScanSum(const double rLocalValue) const
    {
        return rLocalValue;
    }

    virtual void ScanSum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rPartialSums) const
    {}

    virtual void ScanSum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rPartialSums) const
    {}

    // Sendrecv operations

    virtual void SendRecv(
        const std::vector<int>& rSendValues, const int SendDestination,
        std::vector<int>& rRecvValues, const int RecvSource) const
    {}

    virtual void SendRecv(
        const std::vector<double>& rSendValues, const int SendDestination,
        std::vector<double>& rRecvValues, const int RecvSource) const
    {}

    // Broadcast

    virtual void Broadcast(
        int& rBuffer,
        const int SourceRank) const
    {}

    virtual void Broadcast(
        double& rBuffer,
        const int SourceRank) const
    {}

    virtual void Broadcast(
        std::vector<int>& rBuffer,
        const int SourceRank) const
    {}

    virtual void Broadcast(
        std::vector<double>& rBuffer,
        const int SourceRank) const
    {}

    // Scatter operations

    virtual void Scatter(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const int SourceRank) const
    {}

    virtual void Scatter(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const int SourceRank) const
    {}

    virtual void Scatterv(
        const std::vector<int>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<int>& rRecvValues,
        const int SourceRank) const
    {}

    virtual void Scatterv(
        const std::vector<double>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<double>& rRecvValues,
        const int SourceRank) const
    {}

    // Gather operations

    /// Wrapper for MPI_Gather calls (int version).
    /** @param rSendValues Values to be gathered from this rank.
     *  @param rRecvValues Container for the result of the MPI_Allgather call.
     *  @param DestinationRank The rank where the values will be gathered.
     *  @note rRecvValues is only meaningful on rank DestinationRank.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void Gather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const int DestinationRank) const
    {}

    /// Wrapper for MPI_Gather calls (double version).
    /** @param rSendValues Values to be gathered from this rank.
     *  @param rRecvValues Container for the result of the MPI_Allgather call.
     *  @param DestinationRank The rank where the values will be gathered.
     *  @note rRecvValues is only meaningful on rank DestinationRank.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void Gather(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const int DestinationRank) const
    {}

    virtual void Gatherv(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const int DestinationRank) const
    {}

    virtual void Gatherv(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const int DestinationRank) const
    {}

    /// Wrapper for MPI_Allgather calls (int version).
    /** @param rSendValues Values to be gathered from this rank.
     *  @param rRecvValues Container for the result of the MPI_Allgather call.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void AllGather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues) const
    {}

    /// Wrapper for MPI_Allgather calls (double version).
    /** @param rSendValues Values to be gathered from this rank.
     *  @param rRecvValues Container for the result of the MPI_Allgather call.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void AllGather(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues) const
    {}

    ///@}
    ///@name Inquiry
    ///@{

    /// Retrun the parallel rank for this DataCommunicator.
    /** This is a wrapper for calls to MPI_Comm_rank. */
    virtual int Rank() const
    {
        return 0;
    }

    /// Retrun the parallel size for this DataCommunicator.
    /** This is a wrapper for calls to MPI_Comm_size. */
    virtual int Size() const
    {
        return 1;
    }

    /// Check whether this DataCommunicator is aware of parallelism.
    virtual bool IsDistributed() const
    {
        return false;
    }

    ///@}
    ///@name Helper functions for error checking in MPI
    ///@{

    /// This function throws an error on ranks != Sourcerank if Condition evaluates to true.
    /** This method is intended as a helper function to force processes to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF( data_communicator_instance.BroadcastErrorIfTrue(Condition, Root) )
     *  << "Detailed error message in Root rank";
     *
     *  If an error is detected, processes other than Root will fail with a generic error message.
     *  Failing on the Root rank is left to the caller, so that a detailed error message can be
     *  produced.
     *
     *  @note: This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @param SourceRank The rank where the condition is meaningful.
     *  @return The result of evaluating Condition.
     */
    virtual bool BroadcastErrorIfTrue(bool Condition, const int SourceRank) const
    {
        return Condition;
    }

    /// This function throws an error on ranks != Sourcerank if Condition evaluates to false.
    /** This method is intended as a helper function to force processes to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF_NOT( data_communicator_instance.BroadcastErrorIfFalse(Condition, Root) )
     *  << "Detailed error message in Root rank";
     *
     *  If an error is detected, processes other than Root will fail with a generic error message.
     *  Failing on the Root rank is left to the caller, so that a detailed error message can be
     *  produced.
     *
     *  @note: This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @param SourceRank The rank where the condition is meaningful.
     *  @return The result of evaluating Condition.
     */
    virtual bool BroadcastErrorIfFalse(bool Condition, const int SourceRank) const
    {
        return Condition;
    }

    /// This function throws an error on ranks where Condition evaluates to false, if it evaluated to true on a different rank.
    /** This method is intended as a helper function to force processes to stop after an error
     *  is detected on one or more ranks.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF( data_communicator_instance.ErrorIfTrueOnAnyRank(Condition) )
     *  << "Detailed error message in ranks where Condition == true.";
     *
     *  If an error is detected, processes other than those where it was detected will fail with
     *  a generic error message.
     *  Failing on the ranks where the condition is true is left to the caller,
     *  so that a detailed error message can be produced.
     *
     *  @note: This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @return The result of evaluating Condition.
     */
    virtual bool ErrorIfTrueOnAnyRank(bool Condition) const
    {
        return Condition;
    }

    /// This function throws an error on ranks where Condition evaluates to true, if it evaluated to false on a different rank.
    /** This method is intended as a helper function to force processes to stop after an error
     *  is detected on one or more ranks.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF_NOT( data_communicator_instance.ErrorIfFalseOnAnyRank(Condition) )
     *  << "Detailed error message in ranks where Condition == false.";
     *
     *  If an error is detected, processes other than those where it was detected will fail with
     *  a generic error message.
     *  Failing on the ranks where the condition is false is left to the caller,
     *  so that a detailed error message can be produced.
     *
     *  @note: This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @return The result of evaluating Condition.
     */
    virtual bool ErrorIfFalseOnAnyRank(bool Condition) const
    {
        return Condition;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        PrintInfo(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DataCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream
        << "Serial do-nothing version of the Kratos wrapper for MPI communication.\n"
        << "Rank 0 of 1 assumed." << std::endl;
    }

    ///@}

  private:

    ///@name Un accessible methods
    ///@{

    /// Copy constructor.
    DataCommunicator(DataCommunicator const &rOther) = delete;

    /// Assignment operator.
    DataCommunicator &operator=(DataCommunicator const &rOther) = delete;

    ///@}

}; // Class DataCommunicator

template class KRATOS_API(KRATOS_CORE) KratosComponents<DataCommunicator >;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, DataCommunicator const& ThisComponent);

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DataCommunicator &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DataCommunicator &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_DATA_COMMUNICATOR_H_INCLUDED  defined
