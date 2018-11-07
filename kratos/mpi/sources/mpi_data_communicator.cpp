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

#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos {

namespace Internals {

// MPI_Datatype wrappers

template<> inline MPI_Datatype MPIDatatype<int>(const int&)
{
    return MPI_INT;
}

template<> inline MPI_Datatype MPIDatatype<double>(const double&)
{
    return MPI_DOUBLE;
}

template<> inline MPI_Datatype MPIDatatype<array_1d<double,3>>(const array_1d<double,3>&)
{
    return MPI_DOUBLE;
}

}

// MPIDataCommunicator implementation

// Reduce operations

void MPIDataCommunicator::Sum(int& rValue, const int Root) const
{
    int local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, Root, mComm);
}

void MPIDataCommunicator::Sum(double& rValue, const int Root) const
{
    double local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, Root, mComm);
}

void MPIDataCommunicator::Sum(array_1d<double,3>& rValue, const int Root) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Reduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_SUM, Root, mComm);
}

void MPIDataCommunicator::Min(int& rValue, const int Root) const
{
    int local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, Root, mComm);
}

void MPIDataCommunicator::Min(double& rValue, const int Root) const
{
    double local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, Root, mComm);
}

void MPIDataCommunicator::Min(array_1d<double,3>& rValue, const int Root) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Reduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MIN, Root, mComm);
}

void MPIDataCommunicator::Max(int& rValue, const int Root) const
{
    int local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, Root, mComm);
}

void MPIDataCommunicator::Max(double& rValue, const int Root) const
{
    double local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, Root, mComm);
}

void MPIDataCommunicator::Max(array_1d<double,3>& rValue, const int Root) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Reduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MAX, Root, mComm);
}

// Allreduce operations

void MPIDataCommunicator::SumAll(int& rValue) const
{
    int local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, mComm);
}

void MPIDataCommunicator::SumAll(double& rValue) const
{
    double local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, mComm);
}

void MPIDataCommunicator::SumAll(array_1d<double,3>& rValue) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Allreduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_SUM, mComm);
}

void MPIDataCommunicator::MinAll(int& rValue) const
{
    int local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, mComm);
}

void MPIDataCommunicator::MinAll(double& rValue) const
{
    double local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, mComm);
}

void MPIDataCommunicator::MinAll(array_1d<double,3>& rValue) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Allreduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MIN, mComm);
}

void MPIDataCommunicator::MaxAll(int& rValue) const
{
    int local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, mComm);
}

void MPIDataCommunicator::MaxAll(double& rValue) const
{
    double local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, mComm);
}

void MPIDataCommunicator::MaxAll(array_1d<double,3>& rValue) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Allreduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MAX, mComm);
}

}