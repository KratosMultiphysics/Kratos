//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Michael Andre
//                   Philipp Bucher
//

// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // required for the automatic conversion std::vector <=> python-list, see below

// Project includes
#include "mpi_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(mpipython, m)
{
    namespace py = pybind11;

    py::class_<PythonMPIComm>(m,"PythonMPIComm")
    .def(py::init<>())
    .def("barrier",&PythonMPIComm::barrier)
    ;

    typedef int(PythonMPI::*RankFuncType)();
    typedef int(PythonMPI::*SizeFuncType)();
    RankFuncType FRank = &PythonMPI::rank;
    SizeFuncType FSize = &PythonMPI::size;

    // note that for the functions returning a vector the conversion to a python-list is automatically
    // done by pybind, see https://github.com/pybind/pybind11/blob/master/docs/advanced/cast/stl.rst

    const auto py_mpi = py::class_<PythonMPI>(m,"PythonMPI")
    .def_property_readonly("rank",FRank)
    .def_property_readonly("size",FSize)
    .def_property_readonly("world",&PythonMPI::GetWorld,py::return_value_policy::reference_internal )

    .def("broadcast_double",&PythonMPI::broadcast<double>)
    .def("broadcast_int",&PythonMPI::broadcast<int>)

    .def("max_double",&PythonMPI::max<double>)
    .def("max_int",&PythonMPI::max<int>)
    .def("min_double",&PythonMPI::min<double>)
    .def("min_int",&PythonMPI::min<int>)
    .def("sum_double",&PythonMPI::sum<double>)
    .def("sum_int",&PythonMPI::sum<int>)

    .def("max_all_double",&PythonMPI::max_all<double>)
    .def("max_all_int",&PythonMPI::max_all<int>)
    .def("min_all_double",&PythonMPI::min_all<double>)
    .def("min_all_int",&PythonMPI::min_all<int>)
    .def("sum_all_double",&PythonMPI::sum_all<double>)
    .def("sum_all_int",&PythonMPI::sum_all<int>)

    .def("scatter_double", &PythonMPI::scatter<double>)
    .def("scatter_int", &PythonMPI::scatter<int>)

    .def("scatterv_double", &PythonMPI::scatterv<double>)
    .def("scatterv_int", &PythonMPI::scatterv<int>)

    .def("gather_double", &PythonMPI::gather<double>)
    .def("gather_int", &PythonMPI::gather<int>)

    .def("gatherv_double", &PythonMPI::gatherv<double>)
    .def("gatherv_int", &PythonMPI::gatherv<int>)

    .def("allgather_double",&PythonMPI::allgather<double>)
    .def("allgather_int",&PythonMPI::allgather<int>)
    ;

    m.def("GetMPIInterface",&GetMPIInterface,py::return_value_policy::reference);
}

} // Namespace Python
} // Namespace Kratos
