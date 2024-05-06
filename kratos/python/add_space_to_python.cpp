//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_space_to_python.h"
#include "spaces/ublas_space.h"

namespace Kratos::Python
{
void AddSpaceToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;

    py::class_<SparseSpaceType, SparseSpaceType::Pointer>(m,"SparseSpace")
    .def(py::init<>())
    .def_static("Size", &SparseSpaceType::Size)
    .def_static("Size1", &SparseSpaceType::Size1)
    .def_static("Size2", &SparseSpaceType::Size2)
    .def_static("IsDistributed", &SparseSpaceType::IsDistributed)
    .def_static("FastestDirectSolverList", &SparseSpaceType::FastestDirectSolverList)
    ;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    py::class_<LocalSpaceType, LocalSpaceType::Pointer>(m,"LocalSpace")
    .def(py::init<>())
    .def_static("Size", &LocalSpaceType::Size)
    .def_static("Size1", &LocalSpaceType::Size1)
    .def_static("Size2", &LocalSpaceType::Size2)
    .def_static("IsDistributed", &LocalSpaceType::IsDistributed)
    .def_static("FastestDirectSolverList", &LocalSpaceType::FastestDirectSolverList)
    ;
}
}  // namespace Kratos::Python.