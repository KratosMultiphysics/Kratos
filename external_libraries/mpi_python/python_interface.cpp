
#include "mpi_python.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Kratos
{

namespace Python
{

using namespace pybind11;

std::vector<int> (PythonMPI::*gather_int)(PythonMPIComm&, const int, const int) = &PythonMPI::gather;
std::vector<double> (PythonMPI::*gather_double)(PythonMPIComm&, const double, const int) = &PythonMPI::gather;
std::vector<std::vector<int>> (PythonMPI::*gather_list_int)(PythonMPIComm&, const std::vector<int>&, const int) = &PythonMPI::gather;
std::vector<std::vector<double>> (PythonMPI::*gather_list_double)(PythonMPIComm&, const std::vector<double>&, const int) = &PythonMPI::gather;

PYBIND11_MODULE(mpipython, m)
{
    class_<PythonMPIComm>(m,"PythonMPIComm")
    .def(init<>())
    .def("barrier",&PythonMPIComm::barrier)
    ;

    typedef int(PythonMPI::*RankFuncType)();
    typedef int(PythonMPI::*SizeFuncType)();
    RankFuncType FRank = &PythonMPI::rank;
    SizeFuncType FSize = &PythonMPI::size;

    class_<PythonMPI>(m,"PythonMPI")
    .def_property_readonly("rank",FRank)
    .def_property_readonly("size",FSize)
    .def("gather", gather_int)
    .def("gather", gather_double)
    .def("gather", gather_list_int)
    .def("gather", gather_list_double)
    .def("allgather",&PythonMPI::allgather<double>)
    .def("allgather",&PythonMPI::allgather<int>)
    .def_property_readonly("world",&PythonMPI::GetWorld,return_value_policy::reference_internal )
    ;

    m.def("GetMPIInterface",&GetMPIInterface,return_value_policy::reference);
}

} // Namespace Python

} // Namespace Kratos
