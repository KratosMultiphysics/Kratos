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

    /// Copy constructor.
    DataCommunicator(DataCommunicator const &rOther) {};

    /// Destructor.
    virtual ~DataCommunicator() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Barrier() const {}

    // Reduce operations

    virtual int Sum(int& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual double Sum(double& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> Sum(array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual int Min(int& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual double Min(double& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> Min(array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual int Max(int& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual double Max(double& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> Max(array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    // Allreduce operations

    virtual int SumAll(int& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double SumAll(double& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> SumAll(array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual int MinAll(int& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double MinAll(double& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> MinAll(array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual int MaxAll(int& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double MaxAll(double& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual array_1d<double,3> MaxAll(array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    // Scan operations

    virtual int ScanSum(const int& rLocalValue) const
    {
        return rLocalValue;
    }

    virtual double ScanSum(const double& rLocalValue) const
    {
        return rLocalValue;
    }

    // Sendrecv operations

    virtual void SendRecv(
        const std::vector<int>& rSendValues, const unsigned int SendDestination,
        std::vector<int>& rRecvValues, const unsigned int RecvSource) const
    {}

    virtual void SendRecv(
        const std::vector<double>& rSendValues, const unsigned int SendDestination,
        std::vector<double>& rRecvValues, const unsigned int RecvSource) const
    {}

    // Scatter operations

    virtual void Scatter(
        const std::vector<int>& rSendValues,
        std::vector<int>& rRecvValues,
        const unsigned int SourceRank) const
    {}

    virtual void Scatter(
        const std::vector<double>& rSendValues,
        std::vector<double>& rRecvValues,
        const unsigned int SourceRank) const
    {}

    ///@}
    ///@name Inquiry
    ///@{

    virtual int Rank() const
    {
        return 0;
    }

    virtual int Size() const
    {
        return 1;
    }

    virtual bool IsDistributed() const
    {
        return false;
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
    virtual void PrintData(std::ostream &rOStream) const {}

    ///@}

  private:

    ///@name Un accessible methods
    ///@{

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
