//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "add_data_communicator_to_python.h"
#include "includes/define_python.h"
#include "includes/data_communicator.h"

namespace Kratos {

namespace Python {

template<class TValue>
std::vector<TValue> VectorReduceWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pMethod)(const std::vector<TValue>&, std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rLocalValues,
    const int Root)
{
    std::vector<TValue> reduced_values;
    if (rSelf.Rank() == Root)
    {
        reduced_values.resize(rLocalValues.size());
    }
    (rSelf.*pMethod)(rLocalValues, reduced_values, Root);
    return reduced_values;
}

// This wrapper covers all operations that take a std::vector and output to a different std::vector of the same size.
template<class TValue>
std::vector<TValue> VectorBufferTransferWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pMethod)(const std::vector<TValue>&, std::vector<TValue>&) const,
    const std::vector<TValue>& rLocalValues)
{
    std::vector<TValue> reduced_values(rLocalValues.size());
    (rSelf.*pMethod)(rLocalValues, reduced_values);
    return reduced_values;
}

template<class TValue>
std::vector<TValue> VectorSendRecvWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pMethod)(const std::vector<TValue>&, const int, std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rSendValues,
    const int SendDestination,
    const int RecvSource)
{
    std::vector<int> send_size(1, rSendValues.size());
    std::vector<int> recv_size{0};
    rSelf.SendRecv(send_size, SendDestination, recv_size, RecvSource);

    std::vector<TValue> recv_values(recv_size[0]);
    (rSelf.*pMethod)(rSendValues, SendDestination, recv_values, RecvSource);
    return recv_values;
}

template<class TValue>
std::vector<TValue> VectorBroadcastWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pBroadcastMethod)(std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rSourceValues,
    const int SourceRank)
{
    int rank = rSelf.Rank();
    int message_size = rSourceValues.size();
    rSelf.Broadcast(message_size,SourceRank);

    std::vector<TValue> buffer(message_size);
    if (rank == SourceRank)
    {
        buffer = rSourceValues;
    }

    (rSelf.*pBroadcastMethod)(buffer, SourceRank);
    return buffer;
}

template<class TValue>
std::vector<TValue> VectorScatterWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pScatterMethod)(const std::vector<TValue>&, std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rSourceValues,
    const int SourceRank)
{
    const int send_size = rSourceValues.size();
    const int world_size = rSelf.Size();
    KRATOS_ERROR_IF_NOT( send_size % world_size == 0 )
    << "Error in DataCommunicator.Scatter: A message of size " << send_size
    << " cannot be evenly distributed amongst " << world_size << " ranks." << std::endl;
    int message_size = send_size / world_size;

    rSelf.Broadcast(message_size, SourceRank);

    std::vector<TValue> message(message_size);
    (rSelf.*pScatterMethod)(rSourceValues,message,SourceRank);
    return message;
}

template<class TValue>
std::vector<TValue> VectorScattervWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pScattervMethod)(const std::vector<TValue>&, const std::vector<int>&, const std::vector<int>&, std::vector<TValue>&,const int) const,
    const std::vector<std::vector<TValue>>& rSendValues,
    const int SourceRank)
{
    std::vector<TValue> message;
    std::vector<int> message_lenghts;
    std::vector<int> message_offsets;

    if (rSelf.IsDistributed() && (rSelf.Rank() == SourceRank) )
    {
        unsigned int size = rSelf.Size();
        KRATOS_ERROR_IF_NOT(rSendValues.size() == size)
        << "Error in DataCommunicator.Scatterv: expected " << size << " vectors as input, got " << rSendValues.size() << "." << std::endl;

        message_lenghts.resize(size);
        message_offsets.resize(size);
        unsigned int message_size = 0;
        for (unsigned int i = 0; i < rSendValues.size(); i++)
        {
            message_offsets[i] = message_size;
            unsigned int rank_size = rSendValues[i].size();
            message_lenghts[i] = rank_size;
            message_size += rank_size;
        }

        message.resize(message_size);
        for (unsigned int i = 0, counter = 0; i < rSendValues.size(); i++)
        {
            for (unsigned int j = 0; j < rSendValues[i].size(); j++, counter++)
            {
                message[counter] = rSendValues[i][j];
            }
        }
    }

    std::vector<int> recv_size{0};
    rSelf.Scatter(message_lenghts, recv_size, SourceRank);

    std::vector<TValue> recv_message(recv_size[0]);
    rSelf.Scatterv(message, message_lenghts, message_offsets, recv_message, SourceRank);
    return recv_message;
}

template<class TValue>
std::vector<TValue> VectorGatherWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pGatherMethod)(const std::vector<TValue>&, std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rSourceValues,
    const int DestinationRank)
{
    int message_size = rSourceValues.size();
    std::vector<TValue> gathered_values;
    if (rSelf.Rank() == DestinationRank)
    {
        gathered_values.resize(message_size*rSelf.Size());
    }
    (rSelf.*pGatherMethod)(rSourceValues, gathered_values, DestinationRank);
    return gathered_values;
}

template<class TValue>
std::vector<std::vector<TValue>> VectorGathervWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pGathervMethod)(const std::vector<TValue>&, std::vector<TValue>&, const std::vector<int>&, const std::vector<int>&, const int) const,
    const std::vector<TValue>& rSourceValues,
    const int DestinationRank)
{
    std::vector<int> message_sizes_send(1, rSourceValues.size());
    std::vector<int> message_lenghts;
    const int rank = rSelf.Rank();
    const int size = rSelf.Size();
    if (rank == DestinationRank)
    {
        message_lenghts.resize(size);
    }
    rSelf.Gather(message_sizes_send, message_lenghts, DestinationRank);

    std::vector<int> message_offsets;
    std::vector<TValue> gathered_message;
    if (rank == DestinationRank)
    {
        message_offsets.resize(size);
        int message_size = 0;
        for (int i = 0; i < size; i++)
        {
            message_offsets[message_size];
            message_size += message_lenghts[size];
        }
        gathered_message.resize(message_size);
    }
    rSelf.Gatherv(rSourceValues, gathered_message, message_lenghts, message_offsets, DestinationRank);

    std::vector<std::vector<TValue>> gathered_values;
    if (rank == DestinationRank)
    {
        for (int i = 0, counter = 0; i < size; i++)
        {
            gathered_values[i].resize(message_lenghts[i]);
            for (int j = 0; j < message_lenghts[i]; j++, counter++)
            {
                gathered_values[i][j] = gathered_message[counter];
            }
        }
    }
    return gathered_values;
}

template<class TValue>
std::vector<TValue> VectorAllGatherWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pMethod)(const std::vector<TValue>&, std::vector<TValue>&) const,
    const std::vector<TValue>& rLocalValues)
{
    std::vector<TValue> reduced_values(rLocalValues.size()*rSelf.Size());
    (rSelf.*pMethod)(rLocalValues, reduced_values);
    return reduced_values;
}

void AddDataCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<DataCommunicator, DataCommunicator::Pointer>(m,"DataCommunicator")
    .def("Barrier", &DataCommunicator::Barrier)
    // Reduce sum
    .def("Sum", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Sum)
    .def("Sum", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Sum)
    .def("Sum", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Sum)
    .def("SumInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues, const int Root) {
        return VectorReduceWrapper<int>(rSelf, &DataCommunicator::Sum, rLocalValues, Root);
    })
    .def("SumDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues, const int Root) {
        return VectorReduceWrapper<double>(rSelf, &DataCommunicator::Sum, rLocalValues, Root);
    })
    // Reduce min
    .def("Min", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Min)
    .def("Min", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Min)
    .def("Min", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Min)
    .def("MinInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues, const int Root) {
        return VectorReduceWrapper<int>(rSelf, &DataCommunicator::Min, rLocalValues, Root);
    })
    .def("MinDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues, const int Root) {
        return VectorReduceWrapper<double>(rSelf, &DataCommunicator::Min, rLocalValues, Root);
    })
    // Reduce max
    .def("Max", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Max)
    .def("Max", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Max)
    .def("Max", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Max)
    .def("MaxInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues, const int Root) {
        return VectorReduceWrapper<int>(rSelf, &DataCommunicator::Max, rLocalValues, Root);
    })
    .def("MaxDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues, const int Root) {
        return VectorReduceWrapper<double>(rSelf, &DataCommunicator::Max, rLocalValues, Root);
    })
    // Allreduce sum
    .def("SumAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::SumAll)
    .def("SumAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::SumAll)
    .def("SumAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::SumAll)
    .def("SumAllInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        return VectorBufferTransferWrapper<int>(rSelf, &DataCommunicator::SumAll, rLocalValues);
    })
    .def("SumAllDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorBufferTransferWrapper<double>(rSelf, &DataCommunicator::SumAll, rLocalValues);
    })
    // Allreduce min
    .def("MinAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MinAll)
    .def("MinAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MinAll)
    .def("MinAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MinAll)
    .def("MinAllInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        return VectorBufferTransferWrapper<int>(rSelf, &DataCommunicator::MinAll, rLocalValues);
    })
    .def("MinAllDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorBufferTransferWrapper<double>(rSelf, &DataCommunicator::MinAll, rLocalValues);
    })
    // Allreduce max
    .def("MaxAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MaxAll)
    .def("MaxAllInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        return VectorBufferTransferWrapper<int>(rSelf, &DataCommunicator::MaxAll, rLocalValues);
    })
    .def("MaxAllDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorBufferTransferWrapper<double>(rSelf, &DataCommunicator::MaxAll, rLocalValues);
    })
    // ScanSum
    .def("ScanSum", (int (DataCommunicator::*)(const int) const) &DataCommunicator::ScanSum)
    .def("ScanSum", (double (DataCommunicator::*)(const double) const) &DataCommunicator::ScanSum)
    .def("ScanSumInts",[](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        return VectorBufferTransferWrapper<int>(rSelf, &DataCommunicator::ScanSum, rLocalValues);
    })
    .def("ScanSumDoubles",[](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorBufferTransferWrapper<double>(rSelf, &DataCommunicator::ScanSum, rLocalValues);
    })
    // SendRecv
    .def("SendRecvInts",[](DataCommunicator& rSelf, const std::vector<int>& rLocalValues, const int SendDestination, const int RecvSource) {
        return VectorSendRecvWrapper<int>(rSelf, &DataCommunicator::SendRecv, rLocalValues, SendDestination, RecvSource);
    })
    .def("SendRecvDoubles",[](DataCommunicator& rSelf, const std::vector<double>& rLocalValues, const int SendDestination, const int RecvSource) {
        return VectorSendRecvWrapper<double>(rSelf, &DataCommunicator::SendRecv, rLocalValues, SendDestination, RecvSource);
    })
    // Broadcast
    .def("Broadcast", [](DataCommunicator& rSelf, int SourceMessage, const int SourceRank){
        rSelf.Broadcast(SourceMessage,SourceRank);
        return SourceMessage;
    })
    .def("Broadcast",[](DataCommunicator& rSelf, double SourceMessage, const int SourceRank){
        rSelf.Broadcast(SourceMessage,SourceRank);
        return SourceMessage;
    })
    .def("BroadcastInts", [](DataCommunicator& rSelf, const std::vector<int>& rSourceMessage, const int SourceRank) {
        return VectorBroadcastWrapper<int>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    })
    .def("BroadcastDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rSourceMessage, const int SourceRank) {
        return VectorBroadcastWrapper<double>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    })
    // Scatter
    .def("ScatterInts", [](DataCommunicator& rSelf, const std::vector<int>& rSourceMessage, const int SourceRank) {
        return VectorScatterWrapper<int>(rSelf, &DataCommunicator::Scatter, rSourceMessage, SourceRank);
    })
    .def("ScatterDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rSourceMessage, const int SourceRank) {
        return VectorScatterWrapper<double>(rSelf, &DataCommunicator::Scatter, rSourceMessage, SourceRank);
    })
    // Scatterv
    .def("ScattervInts", [](DataCommunicator& rSelf, const std::vector<std::vector<int>>& rSourceMessage, const int SourceRank) {
        return VectorScattervWrapper<int>(rSelf, &DataCommunicator::Scatterv, rSourceMessage, SourceRank);
    })
    .def("ScattervDoubles", [](DataCommunicator& rSelf, const std::vector<std::vector<double>>& rSourceMessage, const int SourceRank) {
        return VectorScattervWrapper<double>(rSelf, &DataCommunicator::Scatterv, rSourceMessage, SourceRank);
    })
    // Gather
    .def("GatherInts",[](DataCommunicator& rSelf, const std::vector<int>& rSourceMessage, const int DestinationRank) {
        return VectorGatherWrapper<int>(rSelf, &DataCommunicator::Gather, rSourceMessage, DestinationRank);
    })
    .def("GatherDoubles",[](DataCommunicator& rSelf, const std::vector<double>& rSourceMessage, const int DestinationRank) {
        return VectorGatherWrapper<double>(rSelf, &DataCommunicator::Gather, rSourceMessage, DestinationRank);
    })
    // Gatherv
    .def("GathervInts",[](DataCommunicator& rSelf, const std::vector<int>& rSourceMessage, const int DestinationRank) {
        return VectorGathervWrapper<int>(rSelf, &DataCommunicator::Gatherv, rSourceMessage, DestinationRank);
    })
    .def("GathervDoubles",[](DataCommunicator& rSelf, const std::vector<double>& rSourceMessage, const int DestinationRank) {
        return VectorGathervWrapper<double>(rSelf, &DataCommunicator::Gatherv, rSourceMessage, DestinationRank);
    })
    // AllGather
    .def("AllGatherInts",[](DataCommunicator& rSelf, const std::vector<int>& rSourceMessage) {
        return VectorAllGatherWrapper<int>(rSelf, &DataCommunicator::AllGather, rSourceMessage);
    })
    .def("AllGatherDoubles",[](DataCommunicator& rSelf, const std::vector<double>& rSourceMessage) {
        return VectorAllGatherWrapper<double>(rSelf, &DataCommunicator::AllGather, rSourceMessage);
    })
    .def("Rank", &DataCommunicator::Rank)
    .def("Size", &DataCommunicator::Size)
    .def("IsDistributed", &DataCommunicator::IsDistributed)
    .def("__str__", PrintObject<DataCommunicator>);
}

} // namespace Python.

} // Namespace Kratos
