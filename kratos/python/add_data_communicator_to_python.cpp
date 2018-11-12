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

// Project includes
#include "add_data_communicator_to_python.h"
#include "includes/define_python.h"
#include "includes/data_communicator.h"

namespace Kratos {

namespace Python {

void AddDataCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<DataCommunicator, DataCommunicator::Pointer>(m,"DataCommunicator")
    .def("Barrier", &DataCommunicator::Barrier)
    .def("Sum", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Sum)
    .def("Sum", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Sum)
    .def("Sum", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Sum)
    .def("Min", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Min)
    .def("Min", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Min)
    .def("Min", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Min)
    .def("Max", (int (DataCommunicator::*)(const int, const int) const) &DataCommunicator::Max)
    .def("Max", (double (DataCommunicator::*)(const double, const int) const) &DataCommunicator::Max)
    .def("Max", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&, const int) const) &DataCommunicator::Max)
    .def("SumAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::SumAll)
    .def("SumAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::SumAll)
    .def("SumAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::SumAll)
    .def("MinAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MinAll)
    .def("MinAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MinAll)
    .def("MinAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MinAll)
    .def("MaxAll", (int (DataCommunicator::*)(const int) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (double (DataCommunicator::*)(const double) const) &DataCommunicator::MaxAll)
    .def("MaxAll", (array_1d<double,3> (DataCommunicator::*)(const array_1d<double,3>&) const) &DataCommunicator::MaxAll)
    .def("ScanSum", (int (DataCommunicator::*)(const int) const) &DataCommunicator::ScanSum)
    .def("ScanSum", (double (DataCommunicator::*)(const double) const) &DataCommunicator::ScanSum)
    .def("Rank", &DataCommunicator::Rank)
    .def("Size", &DataCommunicator::Size)
    .def("__str__", PrintObject<DataCommunicator>);
}

} // namespace Python.

} // Namespace Kratos
