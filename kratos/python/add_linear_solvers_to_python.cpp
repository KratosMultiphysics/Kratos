/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: mossaiby $
//   Date:                $Date: 2008-12-22 14:49:02 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h" 
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "includes/dof.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/scaling_solver.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"
//#include "linear_solvers/superlu_solver.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"



namespace Kratos
{
	
namespace Python
{
  void  AddLinearSolversToPython()
  {
typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;
typedef CGSolver<SpaceType,  LocalSpaceType> CGSolverType;
typedef DeflatedCGSolver<SpaceType,  LocalSpaceType> DeflatedCGSolverType;
typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;
typedef ScalingSolver<SpaceType,  LocalSpaceType> ScalingSolverType;
typedef PowerIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> PowerIterationEigenvalueSolverType;

	bool (LinearSolverType::*pointer_to_solve)(LinearSolverType::SparseMatrixType& rA, LinearSolverType::VectorType& rX, LinearSolverType::VectorType& rB) = &LinearSolverType::Solve;

  using namespace boost::python;

  //****************************************************************************************************
  //preconditioners
  //****************************************************************************************************
	typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
	
    class_<PreconditionerType, PreconditionerType::Pointer>("Preconditioner")
		  .def(self_ns::str(self))
		  ;

	typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
	class_<DiagonalPreconditionerType, DiagonalPreconditionerType::Pointer, bases<PreconditionerType> >("DiagonalPreconditioner")
		  .def(self_ns::str(self))
		  ;

	typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;
	class_<ILUPreconditionerType, ILUPreconditionerType::Pointer, bases<PreconditionerType> >("ILUPreconditioner")
		  .def(self_ns::str(self))
		  ;
	
	typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
	class_<ILU0PreconditionerType, ILU0PreconditionerType::Pointer, bases<PreconditionerType> >("ILU0Preconditioner")
		  .def(self_ns::str(self))
		  ;

  //****************************************************************************************************
  //linear solvers
  //****************************************************************************************************
  class_<LinearSolverType, LinearSolverType::Pointer>("LinearSolver")
		  .def("Initialize",&LinearSolverType::Initialize)
		  .def("Solve",pointer_to_solve)
		  //.def("",&LinearSolverType::)
		  .def(self_ns::str(self))
		  ;

  class_<IterativeSolverType, IterativeSolverType::Pointer, bases<LinearSolverType> >("IterativeSolver")
		  //.def("",&LinearSolverType::)
		  .def(self_ns::str(self))
		  ;

  class_<CGSolverType, CGSolverType::Pointer, bases<IterativeSolverType> >("CGSolver")
		  .def(init<double>())
		  .def(init<double, unsigned int>())
		  .def(init<double, unsigned int,  PreconditionerType::Pointer>())
		  //.def("",&LinearSolverType::)
		  .def(self_ns::str(self))
		  ;

  class_<BICGSTABSolverType, BICGSTABSolverType::Pointer, bases<IterativeSolverType> >("BICGSTABSolver")
		  .def(init<double>())
		  .def(init<double, unsigned int>())
		  .def(self_ns::str(self))
		  .def(init<double, unsigned int,  PreconditionerType::Pointer>())
                  .def("SetTolerance",&BICGSTABSolverType::SetTolerance)
		  ;
  
  class_<TFQMRSolverType, TFQMRSolverType::Pointer, bases<IterativeSolverType> >("TFQMRSolver")
		  .def(init<double>())
		  .def(init<double, unsigned int>())
		  .def(self_ns::str(self))
		  .def(init<double, unsigned int,  PreconditionerType::Pointer>())
		  ;  

    class_<ScalingSolverType, ScalingSolverType::Pointer, bases<LinearSolverType> >("ScalingSolver")
		  .def(init<LinearSolverType::Pointer, bool >())
		  ;

          class_<PowerIterationEigenvalueSolverType, PowerIterationEigenvalueSolverType::Pointer, bases<LinearSolverType> >("PowerIterationEigenvalueSolver")
                  .def(init<double, unsigned int, unsigned int, LinearSolverType::Pointer>())
                          ;

	typedef Reorderer<SpaceType,  LocalSpaceType > ReordererType;
	typedef DirectSolver<SpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
	typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;	

 	class_<ReordererType, ReordererType::Pointer >("Reorderer")
		  .def( init< >() )
		  .def(self_ns::str(self))
		  .def( "Initialize",&ReordererType::Initialize)
		  .def( "Reorder",&ReordererType::Reorder)
		  .def( "InverseReorder",&ReordererType::InverseReorder)
		  ;

 	class_<DirectSolverType, DirectSolverType::Pointer, bases<LinearSolverType> >("DirectSolver")
		  .def( init< >() )
		  .def(self_ns::str(self))
		  ;

	class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, bases<DirectSolverType> >("SkylineLUFactorizationSolver")
		  .def(init< >())
		  .def(self_ns::str(self))
		  ;  

  class_<DeflatedCGSolverType, DeflatedCGSolverType::Pointer, bases<IterativeSolverType> >("DeflatedCGSolver")
		  .def(init<double, int>())
		  .def(init<double, unsigned int, int>())
		  .def(init<double, unsigned int,  PreconditionerType::Pointer, int>())
// 		  .def(init<double, unsigned int,  PreconditionerType::Pointer, ModelPart::Pointer>())
		  //.def("",&LinearSolverType::)
		  .def(self_ns::str(self))
		  ;

//	typedef SuperLUSolver<SparseSpaceType, LocalSpaceType> SuperLUSolverType;
//	class_<SuperLUSolverType, bases<DirectSolverType>, boost::noncopyable >
//		( "SuperLUSolver",
//			init<>() )
  }
	
}  // namespace Python.

} // Namespace Kratos

