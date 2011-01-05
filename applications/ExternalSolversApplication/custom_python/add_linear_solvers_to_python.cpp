//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 12:09:16 $
//   Revision:            $Revision: 1.6 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h" 
#include "spaces/ublas_space.h"
//#include "spaces/parallel_ublas_space.h"

#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "external_includes/superlu_solver.h"
#include "external_includes/superlu_iterative_solver.h"
#include "external_includes/gmres_solver.h"


namespace Kratos
{
    
namespace Python
{
    void  AddLinearSolversToPython()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
	//      typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	// typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;
        typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
        //typedef LinearSolver<ParallelSpaceType,  ParallelLocalSpaceType> ParallelLinearSolverType;
        typedef DirectSolver<SpaceType,  LocalSpaceType> DirectSolverType;
        //typedef Reorderer<ParallelSpaceType,  ParallelLocalSpaceType > ParallelReordererType;
        //typedef DirectSolver<ParallelSpaceType,  ParallelLocalSpaceType, ParallelReordererType > ParallelDirectSolverType;
        typedef SuperLUSolver<SpaceType,  LocalSpaceType> SuperLUSolverType;
	typedef SuperLUIterativeSolver<SpaceType,  LocalSpaceType> SuperLUIterativeSolverType;
        typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
        typedef GMRESSolver<SpaceType, LocalSpaceType> GMRESSolverType;
        typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
        
        using namespace boost::python;
        
        //***************************************************************************
        //linear solvers
        //***************************************************************************
        class_<SuperLUSolverType, bases<DirectSolverType>, boost::noncopyable >
                ( "SuperLUSolver",
                  init<>() );
		
	class_<SuperLUIterativeSolverType, bases<LinearSolverType>, boost::noncopyable >
                ( "SuperLUIterativeSolver",init<>() )
	  .def(init<double,int,int,double,double,double>())
		;
        
        class_<GMRESSolverType, bases<IterativeSolverType>, boost::noncopyable >
                ( "GMRESSolver")
                .def(init<double>())
                .def(init<double, unsigned int>())
                .def(init<double, unsigned int,  PreconditionerType::Pointer>())
                .def(self_ns::str(self))
                ;

  }
	
}  // namespace Python.

} // Namespace Kratos

