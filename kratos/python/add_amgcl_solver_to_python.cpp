//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//


// System includes

// External includes

//this defines are to minimize compilation problems under windows. We could actually use them only when compiling with msvc
#define AMGCL_RUNTIME_DISABLE_MULTICOLOR_GS
#define AMGCL_RUNTIME_DISABLE_PARALLEL_ILU0
#define AMGCL_RUNTIME_DISABLE_SPAI1
#define AMGCL_RUNTIME_DISABLE_CHEBYSHEV

// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "add_amgcl_solver_to_python.h"

#ifndef KRATOS_DISABLE_AMGCL
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"
#endif

namespace Kratos {
namespace Python {

void  AddAMGCLSolverToPython(pybind11::module& m)
{
#ifndef KRATOS_DISABLE_AMGCL
	typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;

    namespace py = pybind11;

    py::enum_<AMGCLSmoother>(m,"AMGCLSmoother")
    .value("SPAI0", SPAI0)
    .value("SPAI1", SPAI1)
    .value("ILU0", ILU0)
    .value("DAMPED_JACOBI",DAMPED_JACOBI)
    .value("GAUSS_SEIDEL",GAUSS_SEIDEL)
    .value("CHEBYSHEV",CHEBYSHEV)
    ;

    py::enum_<AMGCLIterativeSolverType>(m,"AMGCLIterativeSolverType")
    .value("GMRES", GMRES)
    .value("LGMRES", LGMRES)
    .value("FGMRES", FGMRES)
    .value("BICGSTAB", BICGSTAB)
    .value("CG",CG)
    .value("BICGSTAB_WITH_GMRES_FALLBACK",BICGSTAB_WITH_GMRES_FALLBACK)
    .value("BICGSTAB2",BICGSTAB2)
    ;

    py::enum_<AMGCLCoarseningType>(m,"AMGCLCoarseningType")
    .value("RUGE_STUBEN", RUGE_STUBEN)
    .value("AGGREGATION", AGGREGATION)
    .value("SA",SA)
    .value("SA_EMIN",SA_EMIN)
    ;

    typedef AMGCLSolver<SpaceType,  LocalSpaceType> AMGCLSolverType;
    py::class_<AMGCLSolverType,  std::shared_ptr<AMGCLSolverType>, LinearSolverType>
    (m, "AMGCLSolver")
    .def(py::init<AMGCLSmoother,AMGCLIterativeSolverType,double,int,int,int>() )
    .def(py::init<AMGCLSmoother,AMGCLIterativeSolverType,AMGCLCoarseningType ,double,int,int,int, bool>())
    .def(py::init<>())
    .def(py::init<Parameters>())
    .def( "GetResidualNorm",&AMGCLSolverType::GetResidualNorm)
    .def( "GetIterationsNumber",&AMGCLSolverType::GetIterationsNumber)
    ;


    typedef AMGCL_NS_Solver<SpaceType,  LocalSpaceType> AMGCL_NS_SolverType;
    py::class_<AMGCL_NS_SolverType,std::shared_ptr<AMGCL_NS_SolverType>, LinearSolverType >
    (m, "AMGCL_NS_Solver")
    .def(py::init<Parameters>())
    ;
#endif


}

}  // namespace Python.

} // Namespace Kratos

