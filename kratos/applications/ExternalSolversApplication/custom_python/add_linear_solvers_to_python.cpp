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

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "python/add_equation_systems_to_python.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "external_includes/superlu_solver.h"
#include "external_includes/superlu_iterative_solver.h"
#include "external_includes/gmres_solver.h"

#ifdef INCLUDE_FEAST
  #include "external_includes/feast_solver.h"
#endif

#ifndef EXCLUDE_ITSOL
  #include "external_includes/itsol_arms_solver.h"
#endif

#ifdef INCLUDE_PASTIX
  #include "external_includes/pastix_solver.h"
#endif
  

namespace Kratos
{

namespace Python
{
void  AddLinearSolversToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef DirectSolver<SpaceType,  LocalSpaceType> DirectSolverType;
    typedef SuperLUSolver<SpaceType,  LocalSpaceType> SuperLUSolverType;
    typedef SuperLUIterativeSolver<SpaceType,  LocalSpaceType> SuperLUIterativeSolverType;
    typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
    typedef GMRESSolver<SpaceType, LocalSpaceType> GMRESSolverType;
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

    using namespace boost::python;

    //***************************************************************************
    //linear solvers
    //***************************************************************************
#ifdef INCLUDE_FEAST
    typedef FEASTSolver<SpaceType, LocalSpaceType> FEASTSolverType;
    class_<FEASTSolverType, FEASTSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable >
        ( "FEASTSolver", init<Parameters::Pointer>() )
        ;
#endif    
          
    
    class_<SuperLUSolverType, bases<DirectSolverType>, boost::noncopyable >
    ( "SuperLUSolver",
      init<>() )
      .def(init<Parameters>());
      
    class_<SuperLUIterativeSolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "SuperLUIterativeSolver",init<>() )
    .def(init<double,int,int,double,double,double>())
    .def(init<Parameters>())
    ;
    
#ifndef EXCLUDE_ITSOL
    typedef ITSOL_ARMS_Solver<SpaceType,  LocalSpaceType> ITSOL_ARMS_SolverType;
    class_<ITSOL_ARMS_SolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "ITSOL_ARMS_Solver",init<>() )
    .def(init<double,int,int>())
    ;
#endif

#ifdef INCLUDE_PASTIX
    typedef PastixSolver<SpaceType,  LocalSpaceType> PastixSolverType;
    class_<PastixSolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "PastixSolver",init<int,bool>() )
    .def(init<double,int,int,int,bool>())
    .def(init<Parameters>());
    ;
#endif
    
    class_<GMRESSolverType, bases<IterativeSolverType>, boost::noncopyable >
    ( "GMRESSolver")
    .def(init<Parameters >())
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

