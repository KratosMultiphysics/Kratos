//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: pooyan $
//   Date:                $Date: 2008-07-09 13:12:21 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/linear_solver.h"
#include "petsc_application.h"
#include "petsc_solver.h"


namespace Kratos
{
	
namespace Python
{
  void  AddPetscLinearSolversToPython()
  {
	typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
	typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

  	using namespace boost::python;

	typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;
	typedef PetscSolver<SpaceType, LocalSpaceType> PetscSolverType;
	class_<PetscSolverType, bases<LinearSolverType>, boost::noncopyable >
		( "PetscSolver",
			init<>() );
  }
	
}  // namespace Python.

} // Namespace Kratos

