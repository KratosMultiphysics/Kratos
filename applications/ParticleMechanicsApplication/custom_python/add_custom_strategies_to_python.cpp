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
//   Last modified by:    $Author:  ilaria$
//   Date:                $Date:  July 2015$
//   Revision:            $Revision: 1.2 $
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
#include "containers/flags.h"
#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/residual_based_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/MPM_residual_based_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/MPM_strategy.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergence_criteria/displacement_convergence_criterion.hpp"

//schemes
#include "custom_strategies/custom_schemes/residual_based_static_scheme.hpp"
#include "custom_strategies/schemes/MPM_residual_based_bossak_scheme.hpp"


//builders and solvers
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"
//linear solvers
#include "linear_solvers/linear_solver.h"
//#include "structural_application.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

	namespace Python
	{		
		using namespace boost::python;

		void  AddCustomStrategiesToPython()
		{
			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

			//base types
            typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
            typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
            
            //custom strategy types
			typedef MPMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType,2> MPMStrategyType2D;
			typedef MPMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType,3> MPMStrategyType3D;
			
			typedef MPMResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType> MPMResidualBasedNewtonRaphsonStrategyType;
			
			//custom scheme types
			typedef MPMResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >  MPMResidualBasedBossakSchemeType;
			//********************************************************************
			//********************************************************************
// 			class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,	
// 					bases< BaseSolvingStrategyType >,  boost::noncopyable >
// 				("TestStrategy", 
// 				init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
// 				.def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
// 				;
            // MPM Residual Based Bossak Scheme Type
			class_< MPMResidualBasedBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "MPMResidualBasedBossakScheme", init< ModelPart&, double , double >() )

            .def("Initialize", &MPMResidualBasedBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            .def("IterativeExtrapolation", &MPMResidualBasedBossakSchemeType::IterativeExtrapolation)
            ;
            // Strategy Type
			class_< MPMStrategyType2D, bases< BaseSolvingStrategyType >, boost::noncopyable >
				  (  
				  "MPM2D",
				  init < ModelPart&, ModelPart&,  LinearSolverType::Pointer,const Element&, bool , std::string>() )
				
				  
				  .def( "SearchElement", &MPMStrategyType2D::SearchElement)
				  .def( "MP16ShapeFunctions", &MPMStrategyType2D::MP16ShapeFunctions)
				  .def( "MP33ShapeFunctions", &MPMStrategyType2D::MP33ShapeFunctions)
				  ;
			class_< MPMStrategyType3D, bases< BaseSolvingStrategyType >, boost::noncopyable >
				  (  
				  "MPM3D",
				  init < ModelPart&, ModelPart&,  LinearSolverType::Pointer,const Element&, bool, std::string >() )
				
				  
				  .def( "SearchElement", &MPMStrategyType3D::SearchElement)
				  .def( "MP16ShapeFunctions", &MPMStrategyType3D::MP16ShapeFunctions)
				  .def( "MP33ShapeFunctions", &MPMStrategyType3D::MP33ShapeFunctions)
				  ;

			class_< MPMResidualBasedNewtonRaphsonStrategyType,
            bases< BaseSolvingStrategyType >,  boost::noncopyable >
            (
                "MPMResidualBasedNewtonRaphsonStrategy", 
                init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >())

            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
			.def("SetMaxIterationNumber", &MPMResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
			.def("GetMaxIterationNumber", &MPMResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
			.def("SetInitializePerformedFlag", &MPMResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
			.def("GetInitializePerformedFlag", &MPMResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
			.def("SetKeepSystemConstantDuringIterations", &MPMResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
			.def("GetKeepSystemConstantDuringIterations", &MPMResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
			.def("SetFinalizeSolutionStepFlag", &MPMResidualBasedNewtonRaphsonStrategyType::SetFinalizeSolutionStepFlag)
			.def("GetFinalizeSolutionStepFlag", &MPMResidualBasedNewtonRaphsonStrategyType::GetFinalizeSolutionStepFlag)
      ;
    ;


		}

	}  // namespace Python.

} // Namespace Kratos

