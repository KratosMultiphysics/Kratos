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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"


namespace Kratos
{

namespace Python
{
void  AddAMGCLSolverToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;

    using namespace boost::python;

    enum_<AMGCLSmoother>("AMGCLSmoother")
    .value("SPAI0", SPAI0)
    .value("ILU0", ILU0)
    .value("DAMPED_JACOBI",DAMPED_JACOBI)
    .value("GAUSS_SEIDEL",GAUSS_SEIDEL)
    .value("CHEBYSHEV",CHEBYSHEV)
    ;

    enum_<AMGCLIterativeSolverType>("AMGCLIterativeSolverType")
    .value("GMRES", GMRES)
    .value("BICGSTAB", BICGSTAB)
    .value("CG",CG)
    .value("BICGSTAB_WITH_GMRES_FALLBACK",BICGSTAB_WITH_GMRES_FALLBACK)
    .value("BICGSTAB2",BICGSTAB2)
    ;

    enum_<AMGCLCoarseningType>("AMGCLCoarseningType")
    .value("RUGE_STUBEN", RUGE_STUBEN)
    .value("AGGREGATION", AGGREGATION)
    .value("SA",SA)
    .value("SA_EMIN",SA_EMIN)
    ;


    typedef AMGCLSolver<SpaceType,  LocalSpaceType> AMGCLSolverType;
    class_<AMGCLSolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "AMGCLSolver",init<AMGCLSmoother,AMGCLIterativeSolverType,double,int,int,int>() )
    .def(init<AMGCLSmoother,AMGCLIterativeSolverType,AMGCLCoarseningType ,double,int,int,int, bool>())
    .def(init<Parameters>())
    .def( "GetResidualNorm",&AMGCLSolverType::GetResidualNorm)
    .def( "GetIterationsNumber",&AMGCLSolverType::GetIterationsNumber)
    ;


    typedef AMGCL_NS_Solver<SpaceType,  LocalSpaceType> AMGCL_NS_SolverType;
    class_<AMGCL_NS_SolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "AMGCL_NS_Solver", init<Parameters>())
    ;


}

}  // namespace Python.

} // Namespace Kratos

