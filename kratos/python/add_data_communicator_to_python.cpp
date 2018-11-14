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

template<class TValue>
std::vector<TValue> VectorAllReduceWrapper(
    DataCommunicator& rSelf,
    void (DataCommunicator::*pMethod)(const std::vector<TValue>&, std::vector<TValue>&) const,
    const std::vector<TValue>& rLocalValues)
{
    std::vector<TValue> reduced_values(rLocalValues.size());
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
        return VectorAllReduceWrapper<int>(rSelf, &DataCommunicator::SumAll, rLocalValues);
    })
    .def("SumAllDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorAllReduceWrapper<double>(rSelf, &DataCommunicator::SumAll, rLocalValues);
    })
    // Allreduce min
    .def("MinAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MinAll)
    .def("MinAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MinAll)
    .def("MinAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MinAll)
    .def("MinAllInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        return VectorAllReduceWrapper<int>(rSelf, &DataCommunicator::MinAll, rLocalValues);
    })
    .def("MinAllDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorAllReduceWrapper<double>(rSelf, &DataCommunicator::MinAll, rLocalValues);
    })
    // Allreduce max
    .def("MaxAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MaxAll)
    .def("MaxAllInts", [](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        return VectorAllReduceWrapper<int>(rSelf, &DataCommunicator::MaxAll, rLocalValues);
    })
    .def("MaxAllDoubles", [](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        return VectorAllReduceWrapper<double>(rSelf, &DataCommunicator::MaxAll, rLocalValues);
    })
    .def("ScanSum", (int (DataCommunicator::*)(const int) const) &DataCommunicator::ScanSum)
    .def("ScanSum", (double (DataCommunicator::*)(const double) const) &DataCommunicator::ScanSum)
    .def("ScanSumInts",[](DataCommunicator& rSelf, const std::vector<int>& rLocalValues) {
        // I use the same wrapper as in Allreduce, since the two methods have the same signature
        return VectorAllReduceWrapper<int>(rSelf, &DataCommunicator::ScanSum, rLocalValues);
    })
    .def("ScanSumDoubles",[](DataCommunicator& rSelf, const std::vector<double>& rLocalValues) {
        // I use the same wrapper as in Allreduce, since the two methods have the same signature
        return VectorAllReduceWrapper<double>(rSelf, &DataCommunicator::ScanSum, rLocalValues);
    })
    .def("Rank", &DataCommunicator::Rank)
    .def("Size", &DataCommunicator::Size)
    .def("IsDistributed", &DataCommunicator::IsDistributed)
    .def("__str__", PrintObject<DataCommunicator>);
}

} // namespace Python.

} // Namespace Kratos
