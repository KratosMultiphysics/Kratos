
#include "mpi_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

list (PythonMPI::*gather_list)(PythonMPIComm&,list, int) = &PythonMPI::gather;

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
    .def("gather",&PythonMPI::gather<double>)
    .def("gather",&PythonMPI::gather<int>)
    .def("gather", gather_list)
    .def("allgather",&PythonMPI::allgather<double>)
    .def("allgather",&PythonMPI::allgather<int>)
    .def_property_readonly("world",&PythonMPI::GetWorld,return_value_policy::reference_internal )
    ;

    m.def("GetMPIInterface",&GetMPIInterface,return_value_policy::reference);
}

} // Namespace Python

} // Namespace Kratos
