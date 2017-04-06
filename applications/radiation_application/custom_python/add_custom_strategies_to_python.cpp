/*
==============================================================================
KratosConvectionDiffusionApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last modified by:    $Author: julio.marti $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//
//


// System includes 


// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"



//convergence criterias
#include "custom_strategies/strategies/residualbased_convdiffr_strategy.h"
#include "custom_strategies/strategies/residualbased_convdiffr_strategy_nonlinear.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

	namespace Python
	{		
		using namespace boost::python;

		void  AddCustomStrategiesToPython()
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
			typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
			typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

			//********************************************************************
			//********************************************************************
			//


			class_< ResidualBasedConvectionDiffusionrStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedConvectionDiffusionrStrategy", 
				init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int	>() )
				  .def("Clear",&ResidualBasedConvectionDiffusionrStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				;

			class_< ResidualBasedConvectionDiffusionrStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
					bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedConvectionDiffusionrStrategyNonLinear", 
				init<	ModelPart&, LinearSolverType::Pointer,	bool, int, int ,double	>() )
				  .def("Clear",&ResidualBasedConvectionDiffusionrStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >::Clear)
				;	
		}

	}  // namespace Python.

} // Namespace Kratos

