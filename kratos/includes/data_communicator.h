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
#include "containers/array_1d.h"
#include "includes/define.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Serial (do-nothing) version of a wrapper class for MPI communication.
/** @see MPIDataCommunicator for a working distributed memory implementation.
  */
class KRATOS_API(KRATOS_CORE) DataCommunicator
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

    /// Create a new DataCommunicator as a copy of this one.
    /** This method is used in ParallelEnvironment to register DataCommunicators
     *  @see ParallelEnvironment.
     *  @return a unique pointer to the new DataCommunicator.
     */
    virtual DataCommunicator::UniquePointer Clone() const
    {
        return Kratos::make_unique<DataCommunicator>();
    }

    /// Pause program exectution until all threads reach this call.
    /** Wrapper for MPI_Barrier. */
    virtual void Barrier() const {}

    // Reduce operations

    /// Sum rLocalValue across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The summed quantity (meaningful only in Root).
     */
    virtual int Sum(const int rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Sum rLocalValue across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The summed quantity (meaningful only in Root).
     */
    virtual double Sum(const double rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Sum rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The summed quantity (meaningful only in Root).
     */
    virtual array_1d<double,3> Sum(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Sum rLocalValues across all ranks in the Communicator (int vector version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contribution to the sum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The summed quantity (meaningful only in Root).
     */
    virtual std::vector<int> Sum(const std::vector<int>& rLocalValues, const int Root) const
    {
        return rLocalValues;
    }

    /// Sum rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contribution to the sum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The summed quantity (meaningful only in Root).
     */
    virtual std::vector<double> Sum(const std::vector<double>& rLocalValues, const int Root) const
    {
        return rLocalValues;
    }

    /// Sum rLocalValues across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contributions to the sum.
     *  @param[out] rGlobalValues Total sums (meaningful only in Root).
     *  @param[in] Root The rank where the result will be computed.
     */
    virtual void Sum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
    {}

    /// Sum rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contributions to the sum.
     *  @param[out] rGlobalValues Total sums (meaningful only in Root).
     *  @param[in] Root The rank where the result will be computed.
     */
    virtual void Sum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
    {}

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The minimum value (meaningful only in Root).
     */
    virtual int Min(const int rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The minimum value (meaningful only in Root).
     */
    virtual double Min(const double rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The minimum value (meaningful only in Root).
     */
    virtual array_1d<double,3> Min(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the minimum of rLocalValues across all ranks in the Communicator (int vector version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local values to consider in computing the minimum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The minimum values (meaningful only in Root).
     */
    virtual std::vector<int> Min(const std::vector<int>& rLocalValues, const int Root) const
    {
        return rLocalValues;
    }

    /// Obtain the minimum of rLocalValues across all ranks in the Communicator (double vector version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local values to consider in computing the minimum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The minimum values (meaningful only in Root).
     */
    virtual std::vector<double> Min(const std::vector<double>& rLocalValues, const int Root) const
    {
        return rLocalValues;
    }

    /// Obtain the minimum (for each term) of rLocalValues across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contributions to the minimum.
     *  @param[out] rGlobalValues Global minima (meaningful only in Root).
     *  @param[in] Root The rank where the result will be computed.
     */
    virtual void Min(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
    {}

    /// Obtain the minimum (for each term) of rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contributions to the minimum.
     *  @param[out] rGlobalValues Global minima (meaningful only in Root).
     *  @param[in] Root The rank where the result will be computed.
     */
    virtual void Min(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
    {}

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The maximum value (meaningful only in Root).
     */
    virtual int Max(const int rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The maximum value (meaningful only in Root).
     */
    virtual double Max(const double rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The maximum value (meaningful only in Root).
     */
    virtual array_1d<double,3> Max(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValues across all ranks in the Communicator (int vector version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local values to consider in computing the maximum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The maximum values (meaningful only in Root).
     */
    virtual std::vector<int> Max(const std::vector<int>& rLocalValues, const int Root) const
    {
        return rLocalValues;
    }

    /// Obtain the maximum of rLocalValues across all ranks in the Communicator (double vector version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local values to consider in computing the maximum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The maximum values (meaningful only in Root).
     */
    virtual std::vector<double> Max(const std::vector<double>& rLocalValues, const int Root) const
    {
        return rLocalValues;
    }

    /// Obtain the maximum (for each term) of rLocalValues across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contributions to the maximum.
     *  @param[out] rGlobalValues Global maxima (meaningful only in Root).
     *  @param[in] Root The rank where the result will be computed.
     */
    virtual void Max(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
    {}

    /// Obtain the maximum (for each term) of rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValues Local contributions to the maximum.
     *  @param[out] rGlobalValues Global maxima (meaningful only in Root).
     *  @param[in] Root The rank where the result will be computed.
     */
    virtual void Max(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
    {}

    // Allreduce operations

    /// Sum rLocalValue across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Alleduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @return The summed quantity.
     */
    virtual int SumAll(const int rLocalValue) const
    {
        return rLocalValue;
    }

    /// Sum rLocalValue across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Alleduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @return The summed quantity.
     */
    virtual double SumAll(const double rLocalValue) const
    {
        return rLocalValue;
    }

    /// Sum rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Alleduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @return The summed quantity.
     */
    virtual array_1d<double,3> SumAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    /// Sum rLocalValues across all ranks in the Communicator (int vector version).
    /** This is a wrapper to MPI_Alleduce.
     *  @param[in] rLocalValues Local contribution to the sum.
     *  @return The summed quantites.
     */
    virtual std::vector<int> SumAll(const std::vector<int>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Sum rLocalValues across all ranks in the Communicator (double vector version).
    /** This is a wrapper to MPI_Alleduce.
     *  @param[in] rLocalValues Local contribution to the sum.
     *  @return The summed quantities.
     */
    virtual std::vector<double> SumAll(const std::vector<double>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Sum rLocalValues across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local contributions to the sum.
     *  @param[out] rGlobalValues Total sums.
     */
    virtual void SumAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
    {}

    /// Sum rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local contributions to the sum.
     *  @param[out] rGlobalValues Total sums.
     */
    virtual void SumAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
    {}

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @return The minimum value.
     */
    virtual int MinAll(const int rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @return The minimum value.
     */
    virtual double MinAll(const double rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @return The minimum value.
     */
    virtual array_1d<double,3> MinAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the minima of rLocalValues across all ranks in the Communicator (int vector version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local values to consider in computing minima.
     *  @return The minimum values.
     */
    virtual std::vector<int> MinAll(const std::vector<int>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Obtain the minima of rLocalValues across all ranks in the Communicator (double vector version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local values to consider in computing minima.
     *  @return The minimum values.
     */
    virtual std::vector<double> MinAll(const std::vector<double>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Obtain the minimum (for each term) of rLocalValues across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local contributions to the minimum.
     *  @param[out] rGlobalValues Global minima.
     */
    virtual void MinAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
    {}

    /// Obtain the minimum (for each term) of rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local contributions to the minimum.
     *  @param[out] rGlobalValues Global minima.
     */
    virtual void MinAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
    {}

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @return The maximum value.
     */
    virtual int MaxAll(const int rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @return The maximum value.
     */
    virtual double MaxAll(const double rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @return The maximum value.
     */
    virtual array_1d<double,3> MaxAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the maxima of rLocalValues across all ranks in the Communicator (int vector version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local values to consider in computing maxima.
     *  @return The maximum values.
     */
    virtual std::vector<int> MaxAll(const std::vector<int>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Obtain the maxima of rLocalValues across all ranks in the Communicator (double vector version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local values to consider in computing maxima.
     *  @return The maximum values.
     */
    virtual std::vector<double> MaxAll(const std::vector<double>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Obtain the maximum (for each term) of rLocalValues across all ranks in the Communicator (int version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local contributions to the maximum.
     *  @param[out] rGlobalValues Global maxima.
     */
    virtual void MaxAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
    {}

    /// Obtain the maximum (for each term) of rLocalValues across all ranks in the Communicator (double version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValues Local contributions to the maximum.
     *  @param[out] rGlobalValues Global maxima.
     */
    virtual void MaxAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
    {}

    // Scan operations

    /// Compute the partial sums of rLocalValue across all ranks in the Communicator (int version).
    /** The partial sum is the sum of this quantity from rank 0 to the current rank (included).
     *  This is a wrapper to MPI_Scan.
     *  @param[in] rLocalValue Local contribution to the partial sum.
     *  @return The summed quantity.
     */
    virtual int ScanSum(const int rLocalValue) const
    {
        return rLocalValue;
    }

    /// Compute the partial sums of rLocalValue across all ranks in the Communicator (double version).
    /** The partial sum is the sum of this quantity from rank 0 to the current rank (included).
     *  This is a wrapper to MPI_Scan.
     *  @param[in] rLocalValue Local contribution to the partial sum.
     *  @return The summed quantity.
     */
    virtual double ScanSum(const double rLocalValue) const
    {
        return rLocalValue;
    }

    /// Compute the partial sums of rLocalValues across all ranks in the Communicator (int vector version).
    /** The partial sum is the sum of this quantity from rank 0 to the current rank (included).
     *  This is a wrapper to MPI_Scan.
     *  @param[in] rLocalValues Local contributions to the partial sum.
     *  @return The summed quantities.
     */
    virtual std::vector<int> ScanSum(const std::vector<int>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Compute the partial sums of rLocalValues across all ranks in the Communicator (double vector version).
    /** The partial sum is the sum of this quantity from rank 0 to the current rank (included).
     *  This is a wrapper to MPI_Scan.
     *  @param[in] rLocalValues Local contributions to the partial sum.
     *  @return The summed quantities.
     */
    virtual std::vector<double> ScanSum(const std::vector<double>& rLocalValues) const
    {
        return rLocalValues;
    }

    /// Compute the partial sums of rLocalValues across all ranks in the Communicator (int version).
    /** The partial sum is the sum of a quantity from rank 0 to the current rank (included).
     *  This is a wrapper to MPI_Scan.
     *  @param[in] rLocalValues Local contributions to the partial sum.
     *  @param[out] rPartialSums Partial sums for the quantities.
     */
    virtual void ScanSum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rPartialSums) const
    {}

    /// Compute the partial sums of rLocalValues across all ranks in the Communicator (double version).
    /** The partial sum is the sum of a quantity from rank 0 to the current rank (included).
     *  This is a wrapper to MPI_Scan.
     *  @param[in] rLocalValues Local contributions to the partial sum.
     *  @param[out] rPartialSums Partial sums for the quantities.
     */
    virtual void ScanSum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rPartialSums) const
    {}

    // Sendrecv operations

    /// Exchange data with other ranks (int version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[in] RecvSource Rank values are expected from.
     *  @return Received values from rank RecvSource.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<int> SendRecv(
        const std::vector<int>& rSendValues,
        const int SendDestination,
        const int RecvSource) const
    {
        if (SendDestination == RecvSource)
        {
            return rSendValues;
        }
        else
        {
            return std::vector<int>(0);
        }
    }

    /// Exchange data with other ranks (double version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[in] RecvSource Rank values are expected from.
     *  @return Received values from rank RecvSource.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<double> SendRecv(
        const std::vector<double>& rSendValues,
        const int SendDestination,
        const int RecvSource) const
    {
        if (SendDestination == RecvSource)
        {
            return rSendValues;
        }
        else
        {
            return std::vector<double>(0);
        }
    }

    /// Exchange data with other ranks (string version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues String to send to rank SendDestination.
     *  @param[in] SendDestination Rank the string will be sent to.
     *  @param[in] RecvSource Rank the string is expected from.
     *  @return Received string from rank RecvSource.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::string SendRecv(
        const std::string& rSendValues,
        const int SendDestination,
        const int RecvSource) const
    {
        if (SendDestination == RecvSource)
        {
            return rSendValues;
        }
        else
        {
            return std::string("");
        }
    }

    /// Exchange data with other ranks (int version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[out] rRecvValues Received values from rank RecvSource.
     *  @param[in] RecvSource Rank values are expected from.
     *  @param[in] RecvTag Message tag for received values.
     */
    virtual void SendRecv(
        const std::vector<int>& rSendValues, const int SendDestination, const int SendTag,
        std::vector<int>& rRecvValues, const int RecvSource, const int RecvTag) const
    {}

    /// Exchange data with other ranks (double version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[out] rRecvValues Received values from rank RecvSource.
     *  @param[in] RecvSource Rank values are expected from.
     *  @param[in] RecvTag Message tag for received values.
     */
    virtual void SendRecv(
        const std::vector<double>& rSendValues, const int SendDestination, const int SendTag,
        std::vector<double>& rRecvValues, const int RecvSource, const int RecvTag) const
    {}

    /// Exchange data with other ranks (string version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues String to send to rank SendDestination.
     *  @param[in] SendDestination Rank the string will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[out] rRecvValues Received string from rank RecvSource.
     *  @param[in] RecvSource Rank the string is expected from.
     *  @param[in] RecvTag Message tag for received values.
     */
    virtual void SendRecv(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        std::string& rRecvValues, const int RecvSource, const int RecvTag) const
    {}

    // Broadcast

    /// Synchronize a buffer to the value held by the broadcasting rank (int version).
    /** This is a wrapper for MPI_Bcast.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     */
    virtual void Broadcast(
        int& rBuffer,
        const int SourceRank) const
    {}

    /// Synchronize a buffer to the value held by the broadcasting rank (double version).
    /** This is a wrapper for MPI_Bcast.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     */
    virtual void Broadcast(
        double& rBuffer,
        const int SourceRank) const
    {}

    /// Synchronize a buffer to the value held by the broadcasting rank (int version).
    /** This is a wrapper for MPI_Bcast.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     */
    virtual void Broadcast(
        std::vector<int>& rBuffer,
        const int SourceRank) const
    {}

    /// Synchronize a buffer to the value held by the broadcasting rank (double version).
    /** This is a wrapper for MPI_Bcast.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     */
    virtual void Broadcast(
        std::vector<double>& rBuffer,
        const int SourceRank) const
    {}

    // Scatter operations

    /// Wrapper for MPI_Scatter calls (int version).
    /** @param[in] rSendValues Values to be scattered (meaningful only on SourceRank).
     *  @param[in] SourceRank The rank containing the values to be scattered.
     *  @return Scattered values for this rank.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<int> Scatter(
        const std::vector<int>& rSendValues,
        const int SourceRank) const
    {
        if (Rank() == SourceRank)
        {
            return rSendValues;
        }
        else
        {
            return std::vector<int>(0);
        }
    }

    /// Wrapper for MPI_Scatter calls (double version).
    /** @param[in] rSendValues Values to be scattered (meaningful only on SourceRank).
     *  @param[in] SourceRank The rank containing the values to be scattered.
     *  @return Scattered values for this rank.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<double> Scatter(
        const std::vector<double>& rSendValues,
        const int SourceRank) const
    {
        if (Rank() == SourceRank)
        {
            return rSendValues;
        }
        else
        {
            return std::vector<double>(0);
        }
    }

    /// Wrapper for MPI_Scatter calls (int version).
    /** @param[in] rSendValues Values to be scattered (meaningful only on SourceRank).
     *  @param[out] rRecvValues Container for the values to be sent.
     *  @param[in] SourceRank The rank containing the values to be scattered.
     *  @note The expected size of rSendValues is the size of rRecvValues times DataCommunicator::Size().
     */
    virtual void Scatter(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const int SourceRank) const
    {}

    /// Wrapper for MPI_Scatter calls (double version).
    /** @param[in] rSendValues Values to be scattered (meaningful only on SourceRank).
     *  @param[out] rRecvValues Container for the values to be sent.
     *  @param[in] SourceRank The rank containing the values to be scattered.
     *  @note The expected size of rSendValues is the size of rRecvValues times DataCommunicator::Size().
     */
    virtual void Scatter(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const int SourceRank) const
    {}

    // Scatterv operations

    /// Wrapper for MPI_Scatterv calls (int version).
    /** @param[in] rSendValues Values to be scattered (meaningful only on SourceRank).
     *  @param[in] SourceRank The rank containing the values to be scattered.
     *  @return Scattered values for this rank.
     *  @note rSendValues should contain as many vectors as ranks in the communicator.
     *  The i-th vector in the list will be sent to rank i.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<int> Scatterv(
        const std::vector<std::vector<int>>& rSendValues,
        const int SourceRank) const
    {
        int rank = Rank();
        if (rank == SourceRank && rank < static_cast<int>(rSendValues.size()) )
        {
            return rSendValues[rank];
        }
        else
        {
            return std::vector<int>(0);
        }
    }

    /// Wrapper for MPI_Scatterv calls (double version).
    /** @param[in] rSendValues Values to be scattered (meaningful only on SourceRank).
     *  @param[in] SourceRank The rank containing the values to be scattered.
     *  @return Scattered values for this rank.
     *  @note rSendValues should contain as many vectors as ranks in the communicator.
     *  The i-th vector in the list will be sent to rank i.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<double> Scatterv(
        const std::vector<std::vector<double>>& rSendValues,
        const int SourceRank) const
    {
        int rank = Rank();
        if (rank == SourceRank && rank < static_cast<int>(rSendValues.size()) )
        {
            return rSendValues[rank];
        }
        else
        {
            return std::vector<double>(0);
        }
    }

    /// Wrapper for MPI_Scatterv calls (int version).
    /** @param[in] rSendValues Values to be scattered (meaningul only on SourceRank).
     *  @param[in] rSendCounts Number of values to be sent per rank, in order of increasing rank.
     *  @param[in] rSendOffsets Offset from the start of rSendValues of the first value to be sent to each rank.
     *  @param[out] rRecvValues Received values.
     *  The received values at rank i correspond to the range rSendValues[rSendOffsets[i]] to
     *  rSendValues[rSendOffsets[i] + rSendCounts[i]].
     */
    virtual void Scatterv(
        const std::vector<int>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<int>& rRecvValues,
        const int SourceRank) const
    {}

    /// Wrapper for MPI_Scatterv calls (double version).
    /** @param[in] rSendValues Values to be scattered (meaningul only on SourceRank).
     *  @param[in] rSendCounts Number of values to be sent per rank, in order of increasing rank.
     *  @param[in] rSendOffsets Offset from the start of rSendValues of the first value to be sent to each rank.
     *  @param[out] rRecvValues Received values.
     *  The received values at rank i correspond to the range rSendValues[rSendOffsets[i]] to
     *  rSendValues[rSendOffsets[i] + rSendCounts[i]].
     */
    virtual void Scatterv(
        const std::vector<double>& rSendValues,
        const std::vector<int>& rSendCounts,
        const std::vector<int>& rSendOffsets,
        std::vector<double>& rRecvValues,
        const int SourceRank) const
    {}

    // Gather operations

    /// Wrapper for MPI_Gather calls (int version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[in] DestinationRank The rank where the values will be gathered.
     *  @return Gathered values for this rank (meaningful only on DestinationRank).
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<int> Gather(
        const std::vector<int>& rSendValues,
        const int DestinationRank) const
    {
        if (Rank() == DestinationRank)
        {
            return rSendValues;
        }
        else
        {
            return std::vector<int>();
        }
    }

    /// Wrapper for MPI_Gather calls (double version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[in] DestinationRank The rank where the values will be gathered.
     *  @return Gathered values for this rank (meaningful only on DestinationRank).
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<double> Gather(
        const std::vector<double>& rSendValues,
        const int DestinationRank) const
    {
        if (Rank() == DestinationRank)
        {
            return rSendValues;
        }
        else
        {
            return std::vector<double>();
        }
    }

    /// Wrapper for MPI_Gather calls (int version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[out] rRecvValues Container for the result of the MPI_Allgather call.
     *  @param[in] DestinationRank The rank where the values will be gathered.
     *  @note rRecvValues is only meaningful on rank DestinationRank.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void Gather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const int DestinationRank) const
    {}

    /// Wrapper for MPI_Gather calls (double version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[out] rRecvValues Container for the result of the MPI_Allgather call.
     *  @param[in] DestinationRank The rank where the values will be gathered.
     *  @note rRecvValues is only meaningful on rank DestinationRank.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void Gather(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const int DestinationRank) const
    {}

    // Gatherv operations

    /// Wrapper for MPI_Gatherv calls (int version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[in] DestinationRank The rank where the values will be gathered.
     *  @return Gathered values for this rank (meaningful only on DestinationRank).
     *  On DestinationRank, the i-th component of the return corresponds to
     *  the vector received from rank i.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<std::vector<int>> Gatherv(
        const std::vector<int>& rSendValues,
        const int DestinationRank) const
    {
        std::vector<std::vector<int>> output(0);
        if (Rank() == DestinationRank)
        {
            output.resize(1);
            output[0] = rSendValues;
        }
        return output;
    }

    /// Wrapper for MPI_Gatherv calls (double version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[in] DestinationRank The rank where the values will be gathered.
     *  @return Gathered values for this rank (meaningful only on DestinationRank).
     *  On DestinationRank, the i-th component of the return corresponds to
     *  the vector received from rank i.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<std::vector<double>> Gatherv(
        const std::vector<double>& rSendValues,
        const int DestinationRank) const
    {
        std::vector<std::vector<double>> output(0);
        if (Rank() == DestinationRank)
        {
            output.resize(1);
            output[0] = rSendValues;
        }
        return output;
    }

    /// Wrapper for MPI_Gatherv calls (int version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[out] rRecvValues Received values (meaningful only on DestinationRank).
     *  @param[in] rRecvCounts Number of values to be received per rank, in order of increasing rank.
     *  @param[in] rRecvOffsets Offset from the start of rRecvValues of the first value received from each rank.
     *  The gathered values are arranged so that the first rRecvCounts[i] values in rSendValues of rank i
     *  are placed on the range rRecvValues[rRecvOffsets[i]] to rRecvValues[rRecvOffsets[i] + rRecvCounts[i]].
     */
    virtual void Gatherv(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const int DestinationRank) const
    {}

    /// Wrapper for MPI_Gatherv calls (int version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @param[out] rRecvValues Received values (meaningful only on DestinationRank).
     *  @param[in] rRecvCounts Number of values to be received per rank, in order of increasing rank.
     *  @param[in] rRecvOffsets Offset from the start of rRecvValues of the first value received from each rank.
     *  The gathered values are arranged so that the first rRecvCounts[i] values in rSendValues of rank i
     *  are placed on the range rRecvValues[rRecvOffsets[i]] to rRecvValues[rRecvOffsets[i] + rRecvCounts[i]].
     */
    virtual void Gatherv(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const std::vector<int>& rRecvCounts,
        const std::vector<int>& rRecvOffsets,
        const int DestinationRank) const
    {}

    // Allgather operations

    /// Wrapper for MPI_Allgather calls (int version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @return Gathered values.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<int> AllGather(
        const std::vector<int>& rSendValues) const
    {
        return rSendValues;
    }

    /// Wrapper for MPI_Allgather calls (double version).
    /** @param[in] rSendValues Values to be gathered from this rank.
     *  @return Gathered values.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::vector<double> AllGather(
        const std::vector<double>& rSendValues) const
    {
        return rSendValues;
    }

    /// Wrapper for MPI_Allgather calls (int version).
    /** @param rSendValues[in] Values to be gathered from this rank.
     *  @param rRecvValues[out] Container for the result of the MPI_Allgather call.
     *  @note The expected size of rRecvValues is the size of rSendValues times DataCommunicator::Size().
     */
    virtual void AllGather(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues) const
    {}

    /// Wrapper for MPI_Allgather calls (double version).
    /** @param rSendValues[in] Values to be gathered from this rank.
     *  @param rRecvValues[out] Container for the result of the MPI_Allgather call.
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
    ///@name Access
    ///@{

    /// Convenience function to retireve the current default DataCommunicator.
    /** @return A reference to the DataCommunicator instance registered as default in ParallelEnvironment.
     */
    static DataCommunicator& GetDefault();

    ///@}
    ///@name Helper functions for error checking in MPI
    ///@{

    /// This function throws an error on ranks != Sourcerank if Condition evaluates to true.
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF( data_communicator_instance.BroadcastErrorIfTrue(Condition, Root) )
     *  << "Detailed error message in Root rank";
     *
     *  If an error is detected, ranks other than Root will fail with a generic error message.
     *  Failing on the Root rank is left to the caller, so that a detailed error message can be
     *  produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
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
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF_NOT( data_communicator_instance.BroadcastErrorIfFalse(Condition, Root) )
     *  << "Detailed error message in Root rank";
     *
     *  If an error is detected, ranks other than Root will fail with a generic error message.
     *  Failing on the Root rank is left to the caller, so that a detailed error message can be
     *  produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
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
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on one or more ranks.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF( data_communicator_instance.ErrorIfTrueOnAnyRank(Condition) )
     *  << "Detailed error message in ranks where Condition == true.";
     *
     *  If an error is detected, ranks other than those where it was detected will fail with
     *  a generic error message.
     *  Failing on the ranks where the condition is true is left to the caller,
     *  so that a detailed error message can be produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
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
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on one or more ranks.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF_NOT( data_communicator_instance.ErrorIfFalseOnAnyRank(Condition) )
     *  << "Detailed error message in ranks where Condition == false.";
     *
     *  If an error is detected, ranks other than those where it was detected will fail with
     *  a generic error message.
     *  Failing on the ranks where the condition is false is left to the caller,
     *  so that a detailed error message can be produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
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
