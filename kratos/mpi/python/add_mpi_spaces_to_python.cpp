
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:     Philipp Bucher (https://github.com/philbucher)
//


// System includes


// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "mpi/python/add_mpi_spaces_to_python.h"

#include "mpi/spaces/amgcl_mpi_space.h"


namespace Kratos {
namespace Python {

void AddMPISpacesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef amgcl::backend::builtin<double> Backend;
    typedef amgcl::mpi::distributed_matrix<Backend> amgcl_mpi_matrix;
    typedef typename Backend::vector amgcl_mpi_vector;

    typedef AmgclMPISpace<amgcl_mpi_matrix, amgcl_mpi_vector> AmgclMPISparseSpaceType;

    py::class_<AmgclMPISparseSpaceType> (m,"AmgclMPISparseSpace")
        .def(py::init<>())
        ;
}

} // namespace Python.
} // Namespace Kratos
