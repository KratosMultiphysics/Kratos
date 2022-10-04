//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_ginkgo_solvers_to_python.h"
#include "custom_solvers/ginkgo_solver.h"

namespace Kratos {
namespace Python {

void AddGinkgoSolversToPython(pybind11::module& m)
{
	typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;

    namespace py = pybind11;

    typedef GinkgoSolver<SpaceType,  LocalSpaceType> GinkgoSolverType;
    py::class_<GinkgoSolverType, GinkgoSolverType::Pointer, LinearSolverType>(m,"GinkgoSolver")
        .def(py::init< >() )
        .def(py::init<Parameters>())
        .def("__str__", PrintObject<LinearSolverType>)
        ;

}

}  // namespace Python.

} // Namespace Kratos
