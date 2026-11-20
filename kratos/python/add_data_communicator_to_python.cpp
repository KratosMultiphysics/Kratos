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
#include <functional>
#include <memory>

// External includes
#include <pybind11/stl.h>

// Project includes
#include "add_data_communicator_to_python.h"
#include "includes/define_python.h"
#include "includes/data_communicator.h"
#include "includes/parallel_environment.h"

namespace Kratos::Python {

/// Python handle to a non-blocking DataCommunicator operation.
/** Python arguments are converted to temporary C++ objects, which would be destroyed before the operation
 *  completes. This class keeps copies of the communicated buffers alive for the lifetime of the request
 *  and gives access to the received (or reduced) values once the operation is completed.
 */
class PythonDataCommunicatorRequest
{
public:
    PythonDataCommunicatorRequest(
        DataCommunicatorRequest&& rRequest,
        std::shared_ptr<void> pBuffers,
        std::function<pybind11::object()> ResultGetter = nullptr)
        : mpBuffers(std::move(pBuffers)),
          mResultGetter(std::move(ResultGetter)),
          mRequest(std::move(rRequest))
    {
    }

    void Wait()
    {
        mRequest.Wait();
    }

    bool Test()
    {
        return mRequest.Test();
    }

    bool IsCompleted() const
    {
        return mRequest.IsCompleted();
    }

    pybind11::object GetResult() const
    {
        KRATOS_ERROR_IF_NOT(IsCompleted())
        << "The result of a non-blocking operation is only available once it is completed (call Wait or Test first)." << std::endl;
        return mResultGetter ? mResultGetter() : pybind11::none();
    }

private:
    // Declaration order matters: mRequest is destroyed first, so that a pending
    // operation is completed while the buffers it uses are still alive.
    std::shared_ptr<void> mpBuffers;
    std::function<pybind11::object()> mResultGetter;
    DataCommunicatorRequest mRequest;
};

template<class TValue>
PythonDataCommunicatorRequest PythonISend(
    const DataCommunicator& rSelf,
    const TValue& rSendValues,
    const int SendDestination,
    const int SendTag)
{
    auto p_send_values = std::make_shared<TValue>(rSendValues);
    auto request = rSelf.ISend(*p_send_values, SendDestination, SendTag);
    return PythonDataCommunicatorRequest(std::move(request), p_send_values);
}

template<class TValue>
PythonDataCommunicatorRequest PythonIRecv(
    const DataCommunicator& rSelf,
    TValue RecvValues,
    const int RecvSource,
    const int RecvTag)
{
    auto p_recv_values = std::make_shared<TValue>(std::move(RecvValues));
    auto request = rSelf.IRecv(*p_recv_values, RecvSource, RecvTag);
    return PythonDataCommunicatorRequest(std::move(request), p_recv_values, [p_recv_values]() { return pybind11::cast(*p_recv_values); });
}

template<class TValue>
PythonDataCommunicatorRequest PythonIBroadcast(
    const DataCommunicator& rSelf,
    const TValue& rValues,
    const int SourceRank)
{
    auto p_values = std::make_shared<TValue>(rValues);
    auto request = rSelf.IBroadcast(*p_values, SourceRank);
    return PythonDataCommunicatorRequest(std::move(request), p_values, [p_values]() { return pybind11::cast(*p_values); });
}

template<class TValue, class TOperation>
PythonDataCommunicatorRequest PythonIAllReduce(
    const DataCommunicator& rSelf,
    const TValue& rLocalValues,
    TOperation&& rOperation)
{
    // The global buffer is initialized with the local values to have the right size.
    auto p_values = std::make_shared<std::pair<TValue, TValue>>(rLocalValues, rLocalValues);
    auto request = rOperation(rSelf, p_values->first, p_values->second);
    return PythonDataCommunicatorRequest(std::move(request), p_values, [p_values]() { return pybind11::cast(p_values->second); });
}

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
        arg_text = "Int";
    } else if constexpr(std::is_same_v<TDataType, double>) {
        arg_text = "Double";
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
    const std::string& plural_arg_text = (std::is_same_v<TDataType, Matrix> ? "Matrices" : arg_text + "s");
    const std::string& list_of_values = "list_of_" + plural_arg_text;
    const std::string& list_of_v_values = "list_of_" + plural_arg_text + "_per_ranks";

    namespace py = pybind11;

    rDataCommunicatorModule.def("Sum", py::overload_cast<const TDataType&, const int>(&DataCommunicator::Sum, py::const_), py::arg(value_text.c_str()), py::arg("root"));
    rDataCommunicatorModule.def(("Sum" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Sum, py::const_), py::arg(list_of_values.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Min", py::overload_cast<const TDataType&, const int>(&DataCommunicator::Min, py::const_), py::arg(value_text.c_str()), py::arg("root"));
    rDataCommunicatorModule.def(("Min" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Min, py::const_), py::arg(list_of_values.c_str()), py::arg("root"));
    rDataCommunicatorModule.def("Max", py::overload_cast<const TDataType&, const int>(&DataCommunicator::Max, py::const_), py::arg(value_text.c_str()), py::arg("root"));
    rDataCommunicatorModule.def(("Max" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Max, py::const_), py::arg(list_of_values.c_str()), py::arg("root"));

    rDataCommunicatorModule.def("SumAll", py::overload_cast<const TDataType&>(&DataCommunicator::SumAll, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def(("SumAll" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::SumAll, py::const_), py::arg(list_of_values.c_str()));
    rDataCommunicatorModule.def("MinAll", py::overload_cast<const TDataType&>(&DataCommunicator::MinAll, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def(("MinAll" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::MinAll, py::const_), py::arg(list_of_values.c_str()));
    rDataCommunicatorModule.def("MaxAll", py::overload_cast<const TDataType&>(&DataCommunicator::MaxAll, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def(("MaxAll" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::MaxAll, py::const_), py::arg(list_of_values.c_str()));
    if constexpr (std::is_same_v<TDataType, double> || std::is_same_v<TDataType, int>) {
        rDataCommunicatorModule.def("MinLocAll", py::overload_cast<const TDataType&>(&DataCommunicator::MinLocAll, py::const_), py::arg(value_text.c_str()));
        rDataCommunicatorModule.def("MaxLocAll", py::overload_cast<const TDataType&>(&DataCommunicator::MaxLocAll, py::const_), py::arg(value_text.c_str()));
    }

    rDataCommunicatorModule.def("ScanSum", py::overload_cast<const TDataType&>(&DataCommunicator::ScanSum, py::const_), py::arg(value_text.c_str()));
    rDataCommunicatorModule.def(("ScanSum" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::ScanSum, py::const_), py::arg(list_of_values.c_str()));

    rDataCommunicatorModule.def(("SendRecv" + plural_arg_text).c_str(), pybind11::overload_cast<const std::vector<TDataType>&, const int, const int>(&DataCommunicator::SendRecv<std::vector<TDataType>>, pybind11::const_), py::arg(list_of_values.c_str()), py::arg("send_destination"), py::arg("recv_source"));

    rDataCommunicatorModule.def("Broadcast", [](const DataCommunicator& rSelf, TDataType& rSourceMessage, const int SourceRank){
        rSelf.Broadcast(rSourceMessage, SourceRank);
        return rSourceMessage;
    }, py::arg(value_text.c_str()), py::arg("source_rank"));
    rDataCommunicatorModule.def(("Broadcast" + plural_arg_text).c_str(), [](const DataCommunicator& rSelf, std::vector<TDataType>& rSourceMessage, const int SourceRank) {
        return VectorBroadcastWrapper<TDataType>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    }, py::arg(list_of_values.c_str()), py::arg("source_rank"));

    rDataCommunicatorModule.def(("Scatter" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Scatter, py::const_), py::arg(list_of_values.c_str()), py::arg("source_rank"));
    rDataCommunicatorModule.def(("Scatterv" + plural_arg_text).c_str(), py::overload_cast<const std::vector<std::vector<TDataType>>&, const int>(&DataCommunicator::Scatterv, py::const_), py::arg(list_of_v_values.c_str()), py::arg("source_rank"));

    rDataCommunicatorModule.def(("Gather" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Gather, py::const_), py::arg(list_of_values.c_str()), py::arg("destination_rank"));
    rDataCommunicatorModule.def(("Gatherv" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&, const int>(&DataCommunicator::Gatherv, py::const_), py::arg(list_of_values.c_str()), py::arg("destination_rank"));

    rDataCommunicatorModule.def(("AllGather" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::AllGather, py::const_), py::arg(list_of_values.c_str()));
    rDataCommunicatorModule.def(("AllGatherv" + plural_arg_text).c_str(), py::overload_cast<const std::vector<TDataType>&>(&DataCommunicator::AllGatherv, py::const_), py::arg(list_of_values.c_str()));

    rDataCommunicatorModule.def("SynchronizeShape", [](const DataCommunicator& rSelf, TDataType& rValue) { rSelf.SynchronizeShape(rValue); return rValue; }, py::arg(value_text.c_str()));

    // Non-blocking (asynchronous) operations
    if constexpr (std::is_same_v<TDataType, double> || std::is_same_v<TDataType, int>) {
        using VectorType = std::vector<TDataType>;

        rDataCommunicatorModule.def("ISend", &PythonISend<TDataType>, py::arg(value_text.c_str()), py::arg("send_destination"), py::arg("send_tag") = 0);
        rDataCommunicatorModule.def(("ISend" + plural_arg_text).c_str(), &PythonISend<VectorType>, py::arg(list_of_values.c_str()), py::arg("send_destination"), py::arg("send_tag") = 0);

        rDataCommunicatorModule.def(("IRecv" + arg_text).c_str(), [](const DataCommunicator& rSelf, const int RecvSource, const int RecvTag) {
            return PythonIRecv<TDataType>(rSelf, TDataType{}, RecvSource, RecvTag);
        }, py::arg("recv_source"), py::arg("recv_tag") = 0);
        rDataCommunicatorModule.def(("IRecv" + plural_arg_text).c_str(), [](const DataCommunicator& rSelf, const std::size_t NumberOfValues, const int RecvSource, const int RecvTag) {
            return PythonIRecv<VectorType>(rSelf, VectorType(NumberOfValues), RecvSource, RecvTag);
        }, py::arg("number_of_values"), py::arg("recv_source"), py::arg("recv_tag") = 0);

        rDataCommunicatorModule.def("IBroadcast", &PythonIBroadcast<TDataType>, py::arg(value_text.c_str()), py::arg("source_rank"));
        rDataCommunicatorModule.def(("IBroadcast" + plural_arg_text).c_str(), &PythonIBroadcast<VectorType>, py::arg(list_of_values.c_str()), py::arg("source_rank"));

        rDataCommunicatorModule.def("ISumAll", [](const DataCommunicator& rSelf, const TDataType& rLocalValue) {
            return PythonIAllReduce(rSelf, rLocalValue, [](const DataCommunicator& rComm, const TDataType& rLocal, TDataType& rGlobal) { return rComm.ISumAll(rLocal, rGlobal); });
        }, py::arg(value_text.c_str()));
        rDataCommunicatorModule.def(("ISumAll" + plural_arg_text).c_str(), [](const DataCommunicator& rSelf, const VectorType& rLocalValues) {
            return PythonIAllReduce(rSelf, rLocalValues, [](const DataCommunicator& rComm, const VectorType& rLocal, VectorType& rGlobal) { return rComm.ISumAll(rLocal, rGlobal); });
        }, py::arg(list_of_values.c_str()));
        rDataCommunicatorModule.def("IMinAll", [](const DataCommunicator& rSelf, const TDataType& rLocalValue) {
            return PythonIAllReduce(rSelf, rLocalValue, [](const DataCommunicator& rComm, const TDataType& rLocal, TDataType& rGlobal) { return rComm.IMinAll(rLocal, rGlobal); });
        }, py::arg(value_text.c_str()));
        rDataCommunicatorModule.def(("IMinAll" + plural_arg_text).c_str(), [](const DataCommunicator& rSelf, const VectorType& rLocalValues) {
            return PythonIAllReduce(rSelf, rLocalValues, [](const DataCommunicator& rComm, const VectorType& rLocal, VectorType& rGlobal) { return rComm.IMinAll(rLocal, rGlobal); });
        }, py::arg(list_of_values.c_str()));
        rDataCommunicatorModule.def("IMaxAll", [](const DataCommunicator& rSelf, const TDataType& rLocalValue) {
            return PythonIAllReduce(rSelf, rLocalValue, [](const DataCommunicator& rComm, const TDataType& rLocal, TDataType& rGlobal) { return rComm.IMaxAll(rLocal, rGlobal); });
        }, py::arg(value_text.c_str()));
        rDataCommunicatorModule.def(("IMaxAll" + plural_arg_text).c_str(), [](const DataCommunicator& rSelf, const VectorType& rLocalValues) {
            return PythonIAllReduce(rSelf, rLocalValues, [](const DataCommunicator& rComm, const VectorType& rLocal, VectorType& rGlobal) { return rComm.IMaxAll(rLocal, rGlobal); });
        }, py::arg(list_of_values.c_str()));
    }
}


void AddDataCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<PythonDataCommunicatorRequest>(m, "DataCommunicatorRequest")
    .def("Wait", &PythonDataCommunicatorRequest::Wait)
    .def("Test", &PythonDataCommunicatorRequest::Test)
    .def("IsCompleted", &PythonDataCommunicatorRequest::IsCompleted)
    .def("GetResult", &PythonDataCommunicatorRequest::GetResult)
    .def_static("WaitAll", [](const py::list& rRequests) {
        for (const auto& r_request : rRequests) {
            r_request.cast<PythonDataCommunicatorRequest&>().Wait();
        }
    }, py::arg("list_of_requests"))
    ;

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

    data_communicator_module.def("Barrier", &DataCommunicator::Barrier)
    // SendRecv
    .def("SendRecvString",(std::string (DataCommunicator::*)(const std::string&, const int, const int) const) &DataCommunicator::SendRecv)
    // Non-blocking (asynchronous) operations
    .def("IBarrier", [](const DataCommunicator& rSelf) {
        return PythonDataCommunicatorRequest(rSelf.IBarrier(), nullptr);
    })
    .def("ISendString", &PythonISend<std::string>, py::arg("send_string"), py::arg("send_destination"), py::arg("send_tag") = 0)
    .def("IRecvString", [](const DataCommunicator& rSelf, const std::size_t Length, const int RecvSource, const int RecvTag) {
        return PythonIRecv<std::string>(rSelf, std::string(Length, ' '), RecvSource, RecvTag);
    }, py::arg("length"), py::arg("recv_source"), py::arg("recv_tag") = 0)
    // Broadcast
    .def("Broadcast", [](DataCommunicator& rSelf, std::string& rSourceMessage, const int SourceRank){
        rSelf.Broadcast(rSourceMessage, SourceRank);
        return rSourceMessage;
    })
    .def("BroadcastStrings", [](DataCommunicator& rSelf, const std::vector<std::string>& rSourceMessage, const int SourceRank){
        return VectorBroadcastWrapper<std::string>(rSelf, &DataCommunicator::Broadcast, rSourceMessage, SourceRank);
    })
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

} // namespace Kratos::Python.
