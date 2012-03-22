#include <boost/python.hpp>
//#include <boost/python/init.hpp>

#include "mpi_python.h"

namespace Kratos {

namespace Python {

using namespace boost::python;

BOOST_PYTHON_MODULE(mpipython)
{
    class_<PythonMPIComm,boost::noncopyable>("PythonMPIComm")
            .def(init<>())
            .def("barrier",&PythonMPIComm::barrier)
            ;

    typedef int(PythonMPI::*RankFuncType)();
    typedef int(PythonMPI::*SizeFuncType)();
    RankFuncType FRank = &PythonMPI::rank;
    SizeFuncType FSize = &PythonMPI::size;

    class_<PythonMPI,boost::noncopyable>("PythonMPI",no_init) //  init<boost::python::list>())
            .add_property("rank",FRank)
            .add_property("size",FSize)
            .def("gather",&PythonMPI::gather<double>)
            .def("gather",&PythonMPI::gather<int>)
            .def("allgather",&PythonMPI::allgather<double>)
            .def("allgather",&PythonMPI::allgather<int>)
            .add_property("world",make_function(&PythonMPI::GetWorld,return_internal_reference<>() ) )
            ;

    def("GetMPIInterface",&GetMPIInterface,return_value_policy<reference_existing_object>());
}

} // Namespace Python

} // Namespace Kratos
