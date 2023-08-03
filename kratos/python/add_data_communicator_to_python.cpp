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
#include "includes/parallel_environment.h"

namespace Kratos {

namespace Python {

template<class TValue>
std::vector<TValue> VectorBroadcastWrapper(
    const DataCommunicator& rSelf,
    void (DataCommunicator::*pBroadcastMethod)(std::vector<TValue>&, const int) const,
    const std::vector<TValue>& rSourceValues,
    const int SourceRank)
{
    KRATOS_TRY;

    int rank = rSelf.Rank();
    int message_size = rSourceValues.size();
    rSelf.Broadcast(message_size, SourceRank);

    TValue temp{};
    if constexpr(!std::is_same_v<TValue, std::string>) {
        // special case for std::string, since we don't
        // have all the interface support for strings.
        if (rSelf.Rank() == SourceRank && rSourceValues.size() > 0) {
            temp = rSourceValues.front();
        }
        rSelf.SynchronizeShape(temp);
    }

    std::vector<TValue> buffer(message_size, temp);
    if (rank == SourceRank) {
        buffer = rSourceValues;
    }

    (rSelf.*pBroadcastMethod)(buffer, SourceRank);
    return buffer;

    KRATOS_CATCH("")
}

template<class TModuleType, class TDataType>
void AddDataCommunicatorMethodForDataType(
    TModuleType& rDataCommunicatorModule)
{
    std::string arg_text;
    if constexpr(std::is_same_v<TDataType, int>) {
        arg_text = "int";
    } else if constexpr(std::is_same_v<TDataType, double>) {
        arg_text = "double";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        arg_text = "Array3";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 4>>) {
        arg_text = "Array4";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 6>>) {
        arg_text = "Array6";
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 9>>) {
        arg_text = "Array9";
    } else if constexpr(std::is_same_v<TDataType, Vector>) {
        arg_text = "Vector";
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        arg_text = "Matrix";
    } else {
        static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported type.");
    }

    const std::string& value_text = arg_text + "_value";
    const std::string& list_of_values = "list_of_" + arg_text + "_values";
    const std::string& list_of_v_values = "list_of_" + arg_text + "_per_ranks";

    namespace py = pybind11;

    // rDataCommunicatorModule.def("SynchronizeShape", py::overlo)

    rDataCommunicatorModule.def("Sum", py::overload_cast<const TDataType&, const int>(&DataCommunicator::Sum, py::const_), py::arg(value_text.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Sum", py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Sum, py::const_), py::arg(list_of_values.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Min", py::overload_cast<const TDataType&, const int>(&DataCommunicator::Min, py::const_), py::arg(value_text.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Min", py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Min, py::const_), py::arg(list_of_values.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Max", py::overload_cast<const TDataType&, const int>(&DataCommunicator::Max, py::const_), py::arg(value_text.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Max", py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Max, py::const_), py::arg(list_of_values.c_str()), py::arg("root"));

    rDataCommunicatorModule.def("SumAll", py::overload_cast<const TDataType&>(&DataCommunicator::SumAll, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def("SumAll", py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::SumAll, py::const_), py::arg(list_of_values.c_str()));
    rDataCommunicatorModule.def("MinAll", py::overload_cast<const TDataType&>(&DataCommunicator::MinAll, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def("MinAll", py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::MinAll, py::const_), py::arg(list_of_values.c_str()));
    rDataCommunicatorModule.def("MaxAll", py::overload_cast<const TDataType&>(&DataCommunicator::MaxAll, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def("MaxAll", py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::MaxAll, py::const_), py::arg(list_of_values.c_str()));

    rDataCommunicatorModule.def("ScanSum", py::overload_cast<const TDataType&>(&DataCommunicator::ScanSum, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def("ScanSum", py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::ScanSum, py::const_), py::arg(list_of_values.c_str()));

    rDataCommunicatorModule.def("SendRecv", pybind11::overload_cast<const std::vector<TDataType>&, const int, const int>(&DataCommunicator::SendRecv<std::vector<TDataType>>, pybind11::const_), py::arg(list_of_values.c_str()), py::arg("send_destination"), py::arg("recv_source"));

    rDataCommunicatorModule.def("Broadcast", [](const DataCommunicator& rSelf, TDataType& rSourceMessage, const int SourceRank){
        rSelf.Broadcast(rSourceMessage, SourceRank);
        return rSourceMessage;
    }, py::arg(value_text.c_str()), py::arg("source_rank"));

    // TODO: Following vector version of Broadcast cannot be uncommented, because if TDataType=Matrix, then
    //       Pybind11 tries to cast standalone Matrix to std::vector<Matrix> for some reason. Hence, this is deactivated
    //       for the time being.
    // rDataCommunicatorModule.def("Broadcast", [](const DataCommunicator& rSelf, std::vector<TDataType>& rSourceMessage, const int SourceRank) {
    //     return VectorBroadcastWrapper<TDataType>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    // }, py::arg(list_of_values.c_str()), py::arg("source_rank"));

    rDataCommunicatorModule.def("Scatter", py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Scatter, py::const_), py::arg(list_of_values.c_str()), py::arg("source_rank"));
    rDataCommunicatorModule.def("Scatterv", py::overload_cast<const std::vector<std::vector<TDataType>>&, const int>(&DataCommunicator::Scatterv, py::const_), py::arg(list_of_v_values.c_str()), py::arg("source_rank"));

    rDataCommunicatorModule.def("Gather", py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Gather, py::const_), py::arg(list_of_values.c_str()), py::arg("destination_rank"));
    rDataCommunicatorModule.def("Gatherv", py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Gatherv, py::const_), py::arg(list_of_values.c_str()), py::arg("destination_rank"));

    rDataCommunicatorModule.def("AllGather", py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::AllGather, py::const_), py::arg(list_of_values.c_str()));
    rDataCommunicatorModule.def("AllGatherv", py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::AllGatherv, py::const_), py::arg(list_of_values.c_str()));

    rDataCommunicatorModule.def("SynchronizeShape", [](const DataCommunicator& rSelf, TDataType& rValue) { rSelf.SynchronizeShape(rValue); return rValue; }, py::arg(value_text.c_str()));
}


void AddDataCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    auto data_communicator_module = py::class_<DataCommunicator, DataCommunicator::Pointer>(m,"DataCommunicator");

    using data_communicator_module_type = decltype(data_communicator_module);

    AddDataCommunicatorMethodForDataType<data_communicator_module_type, int>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, double>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, array_1d<double, 3>>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, array_1d<double, 4>>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, array_1d<double, 6>>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, array_1d<double, 9>>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, Vector>(data_communicator_module);
    AddDataCommunicatorMethodForDataType<data_communicator_module_type, Matrix>(data_communicator_module);

    // followings  "*Ints" and "*Doubles" versions are kept to have backward compatibility, but they are no longer required.
    data_communicator_module.def("Barrier", &DataCommunicator::Barrier)
    // Reduce sum
    .def("SumInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Sum)
    .def("SumDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Sum)
    // Reduce min
    .def("MinInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Min)
    .def("MinDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Min)
    // Reduce max
    .def("MaxInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int) const) &DataCommunicator::Max)
    .def("MaxDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int) const) &DataCommunicator::Max)
    // Allreduce sum
    .def("SumAllInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::SumAll)
    .def("SumAllDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::SumAll)
    // Allreduce min
    .def("MinAllInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::MinAll)
    .def("MinAllDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::MinAll)
    // Allreduce max
    .def("MaxAllInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::MaxAll)
    .def("MaxAllDoubles", (std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::MaxAll)
    // ScanSum
    .def("ScanSumInts", (std::vector<int> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::ScanSum)
    .def("ScanSumDoubles",(std::vector<double> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::ScanSum)
    // SendRecv
    .def("SendRecvInts",(std::vector<int> (DataCommunicator::*)(const std::vector<int>&, const int, const int) const) &DataCommunicator::SendRecv)
    .def("SendRecvDoubles",(std::vector<double> (DataCommunicator::*)(const std::vector<double>&, const int, const int) const) &DataCommunicator::SendRecv)
    .def("SendRecvString",(std::string (DataCommunicator::*)(const std::string&, const int, const int) const) &DataCommunicator::SendRecv)
    // Broadcast
    .def("Broadcast", [](DataCommunicator& rSelf, std::string& rSourceMessage, const int SourceRank){
        rSelf.Broadcast(rSourceMessage, SourceRank);
        return rSourceMessage;
    })
    .def("BroadcastInts", [](DataCommunicator& rSelf, const std::vector<int>& rSourceMessage, const int SourceRank) {
        return VectorBroadcastWrapper<int>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    })
    .def("BroadcastDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rSourceMessage, const int SourceRank) {
        return VectorBroadcastWrapper<double>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    })
    .def("BroadcastStrings", [](DataCommunicator& rSelf, const std::vector<std::string>& rSourceMessage, const int SourceRank){
        return VectorBroadcastWrapper<std::string>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
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
    // AllGatherv
    .def("AllGathervInts",(std::vector<std::vector<int>> (DataCommunicator::*)(const std::vector<int>&) const) &DataCommunicator::AllGatherv)
    .def("AllGathervDoubles",(std::vector<std::vector<double>> (DataCommunicator::*)(const std::vector<double>&) const) &DataCommunicator::AllGatherv)
    // Common MPI operations
    .def("Rank", &DataCommunicator::Rank)
    .def("Size", &DataCommunicator::Size)
    .def("IsDistributed", &DataCommunicator::IsDistributed)
    .def("IsDefinedOnThisRank", &DataCommunicator::IsDefinedOnThisRank)
    .def("IsNullOnThisRank", &DataCommunicator::IsNullOnThisRank)
    .def_static("GetDefault", []() -> DataCommunicator& {
        KRATOS_WARNING("DataCommunicator") << "This function is deprecated, please retrieve the DataCommunicator through the ModelPart (or by name in special cases)" << std::endl;
        return ParallelEnvironment::GetDefaultDataCommunicator();
    }, py::return_value_policy::reference)
    .def("__str__", PrintObject<DataCommunicator>);
}

} // namespace Python.

} // Namespace Kratos
