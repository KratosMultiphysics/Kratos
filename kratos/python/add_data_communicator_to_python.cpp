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
std::vector<TValue> VectorBroadcastWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pBroadcastMethod)(std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rSourceValues,
    const int SourceRank)
{
    KRATOS_TRY;

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

    KRATOS_CATCH("")
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
    .def("SumInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Sum)
    .def("SumDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Sum)
    // Reduce min
    .def("Min", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Min)
    .def("Min", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Min)
    .def("Min", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Min)
    .def("MinInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Min)
    .def("MinDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Min)
    // Reduce max
    .def("Max", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Max)
    .def("Max", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Max)
    .def("Max", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Max)
    .def("MaxInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Max)
    .def("MaxDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Max)
    // Allreduce sum
    .def("SumAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::SumAll)
    .def("SumAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::SumAll)
    .def("SumAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::SumAll)
    .def("SumAllInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::SumAll)
    .def("SumAllDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::SumAll)
    // Allreduce min
    .def("MinAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MinAll)
    .def("MinAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MinAll)
    .def("MinAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MinAll)
    .def("MinAllInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::MinAll)
    .def("MinAllDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::MinAll)
    // Allreduce max
    .def("MaxAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MaxAll)
    .def("MaxAllInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::MaxAll)
    .def("MaxAllDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::MaxAll)
    // ScanSum
    .def("ScanSum", (int (DataCommunicator::*)(const int) const) &DataCommunicator::ScanSum)
    .def("ScanSum", (double (DataCommunicator::*)(const double) const) &DataCommunicator::ScanSum)
    .def("ScanSumInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::ScanSum)
    .def("ScanSumDoubles",(std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::ScanSum)
    // SendRecv
    .def("SendRecvInts",(std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int, const int) const) &DataCommunicator::SendRecv)
    .def("SendRecvDoubles",(std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int, const int) const) &DataCommunicator::SendRecv)
    .def("SendRecvString",(std::string (DataCommunicator::*)(const std::string&, const int, const int) const) &DataCommunicator::SendRecv)
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
    .def("ScatterInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Scatter)
    .def("ScatterDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Scatter)
    // Scatterv
    .def("ScattervInts", (std::vector<int> (DataCommunicator::*)(const std::vector<std::vector<int>>&, const int) const) &DataCommunicator::Scatterv)
    .def("ScattervDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<std::vector<double>>&, const int) const) &DataCommunicator::Scatterv)
    // Gather
    .def("GatherInts",(std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Gather)
    .def("GatherDoubles",(std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Gather)
    // Gatherv
    .def("GathervInts",(std::vector<std::vector<int>> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Gatherv)
    .def("GathervDoubles",(std::vector<std::vector<double>> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Gatherv)
    // AllGather
    .def("AllGatherInts",(std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::AllGather)
    .def("AllGatherDoubles",(std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::AllGather)
    .def("Rank", &DataCommunicator::Rank)
    .def("Size", &DataCommunicator::Size)
    .def("IsDistributed", &DataCommunicator::IsDistributed)
    .def("IsDefinedOnThisRank", &DataCommunicator::IsDefinedOnThisRank)
    .def("IsNullOnThisRank", &DataCommunicator::IsNullOnThisRank)
    .def_static("GetDefault", &DataCommunicator::GetDefault, py::return_value_policy::reference)
    .def("__str__", PrintObject<DataCommunicator>);
}

} // namespace Python.

} // Namespace Kratos
