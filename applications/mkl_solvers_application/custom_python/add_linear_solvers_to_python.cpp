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
#include "spaces/parallel_ublas_space.h"

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
        typedef ParallelUblasSpace<double, CompressedMatrix, Vector> ParallelSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;
        typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
        typedef LinearSolver<ParallelSpaceType,  ParallelLocalSpaceType> ParallelLinearSolverType;
        typedef DirectSolver<SpaceType,  LocalSpaceType> DirectSolverType;
        typedef Reorderer<ParallelSpaceType,  ParallelLocalSpaceType > ParallelReordererType;
        typedef DirectSolver<ParallelSpaceType,  ParallelLocalSpaceType, ParallelReordererType > ParallelDirectSolverType;
        typedef MKLPardisoSolver<SpaceType, LocalSpaceType> MKLPardisoSolverType;
        typedef MKLPardisoSolver<ParallelSpaceType, ParallelLocalSpaceType> ParallelMKLPardisoSolverType;
        typedef MKLGMRESSolver<SpaceType, LocalSpaceType> MKLGMRESSolverType;
        typedef MKLGMRESSolver<ParallelSpaceType, ParallelLocalSpaceType> ParallelMKLGMRESSolverType;
        typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
        typedef IterativeSolver<ParallelSpaceType, ParallelLocalSpaceType> ParallelIterativeSolverType;
        typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
        
        using namespace boost::python;
        
        //***************************************************************************
        //linear solvers
        //***************************************************************************
        
        class_<MKLPardisoSolverType, MKLPardisoSolverType::Pointer,
        bases<DirectSolverType> >( "MKLPardisoSolver" )
                .def(init<unsigned int>() )
                .def(self_ns::str(self))
                ;
        
        class_<ParallelMKLPardisoSolverType, ParallelMKLPardisoSolverType::Pointer,
        bases<ParallelDirectSolverType> >( "ParallelMKLPardisoSolver" )
                .def(init<unsigned int>() )
                .def(self_ns::str(self))
                ;
        
        class_<MKLGMRESSolverType, MKLGMRESSolverType::Pointer,
        bases<DirectSolverType> >( "MKLGMRESSolver" )
                .def(init<>() )
                .def(self_ns::str(self))
                ;
        
        class_<ParallelMKLGMRESSolverType, ParallelMKLGMRESSolverType::Pointer,
        bases<ParallelDirectSolverType> >( "ParallelMKLGMRESSolver" )
                .def(init<>() )
                .def(self_ns::str(self))
                ;
  }
	
}  // namespace Python.

} // Namespace Kratos

