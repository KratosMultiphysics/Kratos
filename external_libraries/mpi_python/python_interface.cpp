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
    .def("broadcast",&PythonMPI::broadcast<double>)
    .def("broadcast",&PythonMPI::broadcast<int>)
    .def("reduce",&PythonMPI::reduce<double>)
    .def("reduce",&PythonMPI::reduce<int>)
    .def("allreduce",&PythonMPI::allreduce<double>)
    .def("allreduce",&PythonMPI::allreduce<int>)
    .def("gather", &PythonMPI::gather<double>)
    .def("gather", &PythonMPI::gather<int>)
    .def("gatherv", &PythonMPI::gatherv<double>)
    .def("gatherv", &PythonMPI::gatherv<int>)
    .def("allgather",&PythonMPI::allgather<double>)
    .def("allgather",&PythonMPI::allgather<int>)
    .def_property_readonly("world",&PythonMPI::GetWorld,py::return_value_policy::reference_internal )
    ;

    py::enum_<PythonMPI::MPI_Operation>(py_mpi, "MPI_op")
    .value("MAX", PythonMPI::MPI_Operation::MAX)
    .value("MIN", PythonMPI::MPI_Operation::MIN)
    .value("SUM", PythonMPI::MPI_Operation::SUM)
    ;

    m.def("GetMPIInterface",&GetMPIInterface,py::return_value_policy::reference);
}

} // Namespace Python
} // Namespace Kratos
