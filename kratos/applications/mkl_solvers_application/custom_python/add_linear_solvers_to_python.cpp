//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2009-01-14 09:40:45 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h" 
#include "spaces/ublas_space.h"

#ifdef _OPENMP
#include "spaces/parallel_ublas_space.h"
#endif
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "external_includes/mkl_pardiso_solver.h"
#include "external_includes/mkl_gmres_solver.h"


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
        typedef MKLPardisoSolver<SpaceType, LocalSpaceType> MKLPardisoSolverType;
        typedef MKLGMRESSolver<SpaceType, LocalSpaceType> MKLGMRESSolverType;
        typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
        typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

#ifdef _OPENMP
        typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSpaceType;
        typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;
        typedef LinearSolver<ParallelSpaceType,  ParallelLocalSpaceType> ParallelLinearSolverType;
        typedef Reorderer<ParallelSpaceType,  ParallelLocalSpaceType > ParallelReordererType;
        typedef DirectSolver<ParallelSpaceType,  ParallelLocalSpaceType, ParallelReordererType > ParallelDirectSolverType;
        typedef MKLPardisoSolver<ParallelSpaceType, ParallelLocalSpaceType> ParallelMKLPardisoSolverType;
        typedef MKLGMRESSolver<SpaceType, LocalSpaceType> MKLGMRESSolverType;
        typedef MKLGMRESSolver<ParallelSpaceType, ParallelLocalSpaceType> ParallelMKLGMRESSolverType;
        typedef IterativeSolver<ParallelSpaceType, ParallelLocalSpaceType> ParallelIterativeSolverType;
#endif


        using namespace boost::python;

        //***************************************************************************
        //linear solvers
        //***************************************************************************
        
        class_<MKLPardisoSolverType, MKLPardisoSolverType::Pointer,
        bases<DirectSolverType> >( "MKLPardisoSolver" )
                .def(init<unsigned int>() )
                .def(self_ns::str(self))
                ;
        
#ifdef _OPENMP
        class_<ParallelMKLPardisoSolverType, ParallelMKLPardisoSolverType::Pointer,
        bases<ParallelDirectSolverType> >( "ParallelMKLPardisoSolver" )
                .def(init<unsigned int>() )
                .def(self_ns::str(self))
                ;
#endif

        class_<MKLGMRESSolverType, MKLGMRESSolverType::Pointer,
        bases<DirectSolverType> >( "MKLGMRESSolver" )
                .def(init<>() )
                .def(self_ns::str(self))
                ;
        
#ifdef _OPENMP
        class_<ParallelMKLGMRESSolverType, ParallelMKLGMRESSolverType::Pointer,
        bases<ParallelDirectSolverType> >( "ParallelMKLGMRESSolver" )
                .def(init<>() )
                .def(self_ns::str(self))
                ;
#endif
  }
	
}  // namespace Python.

} // Namespace Kratos

