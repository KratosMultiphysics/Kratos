/*
==============================================================================
KratosStructuralApplication 
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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:49:07 $
//   Revision:            $Revision: 1.20 $
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

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/convergencecriterias/multiphaseflow_criteria.h"
//#include "custom_strategies/convergencecriterias/residual_displacement_criteria.h"
//#include "custom_strategies/convergencecriterias/res_dis_criteria.h"
// #include "solving_strategies/convergencecriterias/displacement_criteria.h"
//#include "solving_strategies/convergencecriterias/new_galerkin_displacement_criteria.h"


//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/residualbased_arc_length_strategy.h"
#include "custom_strategies/strategies/residualbased_newton_raphson_line_search_strategy.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
// #include "custom_strategies/strategies/residualbased_newton_raphson_strategy_newmark.h"
// #include "custom_strategies/strategies/residualbased_uzawa_newton_raphson_strategy.h"
// #include "custom_strategies/strategies/residualbased_uzawa_newton_raphson_strategy_newmark.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
// #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/schemes/residualbased_predictorcorrector_bossak_scheme.h"
#include "custom_strategies/schemes/residualbased_central_diferences_scheme.h"
// #include "custom_strategies/schemes/testing_scheme.h"
#include "custom_strategies/schemes/residualbased_predictorcorrector_bossak_scheme_rotation.h"
#include "custom_strategies/schemes/residualbased_predictorcorrector_relaxation_scheme.h"
#include "custom_strategies/schemes/residualbased_newmark_scheme.h"
#include "custom_strategies/schemes/composit_scheme.h"
#include "custom_strategies/schemes/volumetric_scheme.h"

#include "custom_strategies/schemes/inner_volumetric_scheme.h"
#include "custom_strategies/schemes/inner_volumetric_dynamic_scheme.h"
//#include "custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme.h"

//#include "structural_application/custom_strategies/schemes/residualbased_galerkin_scheme.h"

//builder_and_solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
// #include "custom_strategies/builder_and_solvers/modal_analysis_builder_and_solver.h"
//#include "custom_strategies/builder_and_solvers/modal_analysis_builder_and_solver.h"

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
            typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >
                    BaseSolvingStrategyType;
                    
            //typedef ResidualBasedUzawaNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType,
            //LinearSolverType > ResidualBasedUzawaNewtonRaphsonStrategyType;

            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

            typedef ResidualBasedPredictorCorrectorBossakScheme< SparseSpaceType, LocalSpaceType >
                    ResidualBasedPredictorCorrectorBossakSchemeType;
            typedef ResidualBasedPredictorCorrectorBossakRotationScheme< SparseSpaceType, LocalSpaceType >
                    ResidualBasedPredictorCorrectorBossakRotationSchemeType;

            
            typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType;

//             typedef TestingScheme< SparseSpaceType, LocalSpaceType >
//                     TestingSchemeType;
            
            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;
            
            typedef MultiPhaseFlowCriteria< SparseSpaceType,  LocalSpaceType >
                    MultiPhaseFlowCriteriaType;
                    
            typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>
                    BuilderAndSolverType;

//	typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< SparseSpaceType, LocalSpaceType > //ResidualBasedPredictorCorrectorVelocityBossakSchemeType;
                    
//             typedef ModalAnalysisBuilderAndSolver<SparseSpaceType, LocalSpaceType,
//                     LinearSolverType> ModalAnalysisBuilderAndSolverType;
	typedef CompositScheme< SparseSpaceType, LocalSpaceType > CompositSchemeType;
                    
	typedef VolumetricScheme< 2, SparseSpaceType, LocalSpaceType > VolumetricSchemeType2D;
	typedef VolumetricScheme< 3, SparseSpaceType, LocalSpaceType > VolumetricSchemeType3D;

	typedef InnerVolumetricScheme< 2, SparseSpaceType, LocalSpaceType > InnerVolumetricSchemeType2D;
	typedef InnerVolumetricScheme< 3, SparseSpaceType, LocalSpaceType > InnerVolumetricSchemeType3D;

	typedef InnerVolumetricDynamicScheme< 2, SparseSpaceType, LocalSpaceType > InnerVolumetricDynamicSchemeType2D;
	typedef InnerVolumetricDynamicScheme< 3, SparseSpaceType, LocalSpaceType > InnerVolumetricDynamicSchemeType3D;

	
	typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
	//typedef ResidualCriteria < SparseSpaceType, LocalSpaceType >::Pointer TResidual;
	//typedef DisplacementCriteria < SparseSpaceType, LocalSpaceType>::Pointer TDisplacement;

	
// 					;
           //********************************************************************
           //********************************************************************
           class_< ResidualBasedPredictorCorrectorBossakSchemeType,
           bases< BaseSchemeType >,  boost::noncopyable >
                   (
                    "ResidualBasedPredictorCorrectorBossakScheme", init< double >()
                   );

           class_< ResidualBasedCentralDiferencesScheme< SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
                   (
                    "ResidualBasedCentralDiferencesScheme", init< ModelPart&,  double, bool,  bool >())
                   ;

	   class_< ResidualBasedPredictorCorrectorBossakRotationSchemeType,
 	   bases< BaseSchemeType >,  boost::noncopyable >
		    (
		     "ResidualBasedPredictorCorrectorBossakRotationScheme", init< double >()
	            );
    
           typedef ResidualBasedPredictorCorrectorRelaxationScheme< SparseSpaceType,
           LocalSpaceType > ResidualBasedPredictorCorrectorRelaxationSchemeType;
           
           class_< ResidualBasedPredictorCorrectorRelaxationSchemeType,
           bases< BaseSchemeType >,  boost::noncopyable >
                   (
                    "ResidualBasedPredictorCorrectorRelaxationScheme", init< double, double >()
                   );

           class_< ResidualBasedNewmarkSchemeType,
           bases< BaseSchemeType >, boost::noncopyable >
                   (
                    "ResidualBasedNewmarkScheme", init< double >()
                   );

// 			class_< TestingSchemeType,
// 			bases< BaseSchemeType >,  boost::noncopyable >
// 					(
// 					"TestingScheme", init< >()
// 					);
           
           class_< MultiPhaseFlowCriteriaType,
           bases< ConvergenceCriteriaBaseType >, boost::noncopyable >
                   ("MultiPhaseFlowCriteria", init<double, double >() )
                   ;

//            class_< ModalAnalysisBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
//                    (
//                     "ModalAnalysisBuilderAndSolver", init<LinearSolverType::Pointer>()
//                    );



//	class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeType,
//			bases< BaseSchemeType >,  boost::noncopyable >
//					(
//					"ResidualBasedPredictorCorrectorVelocityBossakScheme", init< double >()
//					);
			
			class_< VolumetricSchemeType2D,
			bases< BaseSchemeType >,  boost::noncopyable >
					(
					"VolumetricScheme2D"
					)
			.def("CalculateCauchyStress",&VolumetricSchemeType2D::CalculateCauchyStress)
			;

			class_< InnerVolumetricSchemeType2D,
			bases< BaseSchemeType >,  boost::noncopyable >
					(
					"InnerVolumetricScheme2D"
					);

			class_< InnerVolumetricDynamicSchemeType2D,
			bases< BaseSchemeType >,  boost::noncopyable >
					(
					"InnerVolumetricDynamicScheme2D"
					);

			class_< CompositSchemeType,
			bases< BaseSchemeType >,  boost::noncopyable >
					(
					"CompositScheme", init< BaseSchemeType&, BaseSchemeType& >()
					);
			
			class_< ResidualBasedArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedArcLenghtStrategy", 
				init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, 
                                unsigned int, unsigned int,double,bool, bool, bool,bool
				>() )
				;
			
			
                        
			class_< ResidualBasedNewtonRaphsonLineSearchesStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,bases< BaseSolvingStrategyType >,  boost::noncopyable >
				("ResidualBasedNewtonRaphsonLineSearchesStrategy", 
				init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, unsigned int, unsigned int, double, double, double, double, bool, bool, bool, bool
				>() )
				;




// 			class_<Residual_Displacement_Criteria<SparseSpaceType, LocalSpaceType >,
// 			         bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,  
// 			         boost::noncopyable >
// 			        ("ResidualDisplacementCriteria", init< double, double>() );
// 
// 			class_<ResDisCriteria<SparseSpaceType, LocalSpaceType >,
// 			         bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,  
// 			         boost::noncopyable >
// 			        ("ResDisCriteria", init< TResidual,TDisplacement >());

           /*
           class_< ModalAnalysisBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
                   (
                    "ModalAnalysisBuilderAndSolver", init<LinearSolverType::Pointer>()
                   );
*/
        }
    }  // namespace Python.

} // Namespace Kratos

