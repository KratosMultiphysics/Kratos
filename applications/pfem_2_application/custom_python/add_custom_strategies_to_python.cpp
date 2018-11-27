/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes
//#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
//#include <boost/timer.hpp>

#include <pybind11/pybind11.h>
// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_strategies/explicit_strategy.h"
#include "custom_strategies/pfem_2_monolithic_slip_scheme.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/fracstep_GLS_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

	namespace Python
	{
	  namespace py = pybind11;

		void  AddCustomStrategiesToPython(pybind11::module& m)
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
			typedef ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
			typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

			//********************************************************************
			//********************************************************************
			py::class_< PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType  >::Pointer, BaseSolvingStrategyType >(m,"PFEM2_Explicit_Strategy")
			  .def(py::init<	ModelPart&, int, bool >() )

			  //initialize and finalize. the others are standard and are in the ExplicitStrategy
			  .def("InitializeSolutionStep",&PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeSolutionStep)
			  .def("FinalizeSolutionStep",&PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FinalizeSolutionStep)
			  ;
			py::class_< Fluid_Phase_PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, Fluid_Phase_PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer,BaseSolvingStrategyType >(m,"Fluid_Phase_PFEM2_Explicit_Strategy")
			  .def(py::init<ModelPart&, int, bool >() )

			  //initialize and finalize. the others are standard and are in the ExplicitStrategy
			  .def("InitializeSolutionStep",&Fluid_Phase_PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeSolutionStep)
			  .def("FinalizeSolutionStep",&Fluid_Phase_PFEM2_Explicit_Strategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FinalizeSolutionStep)
			  ;

			py::class_< PFEM2MonolithicSlipScheme< SparseSpaceType, LocalSpaceType >, PFEM2MonolithicSlipScheme< SparseSpaceType, LocalSpaceType>::Pointer>(m,"PFEM2MonolithicSlipScheme")
			  .def(py::init<unsigned int >())// constructor without a turbulence model
			  ;

			py::class_< FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >:: Pointer> (m,"FracStepStrategy")
			  .def(py::init<	ModelPart&, LinearSolverType::Pointer,LinearSolverType::Pointer,bool, double, double,int, int, unsigned int, unsigned int,bool >() )
			  .def("Solve", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Solve)
			  .def("SolveStep1", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep1)
			  .def("SolveStep2", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep2)
			  .def("SolveStep3", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep3)
			  .def("SolveStep4", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep4)
			  .def("IterativeSolve", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::IterativeSolve)
			  .def("SavePressureIteration", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SavePressureIteration)
			  .def("SolvePressure", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolvePressure)
			  .def("SolveStep7", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SolveStep7)
			  .def("FractionalVelocityIteration", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FractionalVelocityIteration)
			  .def("ConvergenceCheck", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ConvergenceCheck)
			  .def("Clear", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
			  .def("Compute", &FracStepStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Compute)
			  ;

		}

	}  // namespace Python.

} // Namespace Kratos


