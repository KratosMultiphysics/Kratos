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
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/convergencecriterias/multiphaseflow_criteria.h"
// #include "custom_strategies/convergencecriterias/resisualbased_multiphase_criteria.h"
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
#include "custom_strategies/strategies/residualbased_central_differences_strategy.h"
//#include "custom_strategies/schemes/residualbased_central_diferences_scheme.h"
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
// #include "custom_strategies/builder_and_solvers/multiphase_builder_and_solver.h"
// #include "custom_strategies/builder_and_solvers/modal_analysis_builder_and_solver.h"
//#include "custom_strategies/builder_and_solvers/modal_analysis_builder_and_solver.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{


namespace Python
{

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;

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

//     typedef ResidualBasedMultiPhaseCriteria< SparseSpaceType, LocalSpaceType > ResidualBasedMultiPhaseCriteriaType;

    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>
    BuilderAndSolverType;

//            typedef MultiPhaseBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> MultiPhaseBuilderAndSolverType;

//	typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< SparseSpaceType, LocalSpaceType > //ResidualBasedPredictorCorrectorVelocityBossakSchemeType;

//             typedef ModalAnalysisBuilderAndSolver<SparseSpaceType, LocalSpaceType,
//                     LinearSolverType> ModalAnalysisBuilderAndSolverType;
    typedef CompositScheme< SparseSpaceType, LocalSpaceType > CompositSchemeType;

    typedef VolumetricScheme< 2, SparseSpaceType, LocalSpaceType > VolumetricSchemeType2D;
    //typedef VolumetricScheme< 3, SparseSpaceType, LocalSpaceType > VolumetricSchemeType3D;

    typedef InnerVolumetricScheme< 2, SparseSpaceType, LocalSpaceType > InnerVolumetricSchemeType2D;
    //typedef InnerVolumetricScheme< 3, SparseSpaceType, LocalSpaceType > InnerVolumetricSchemeType3D;

    typedef InnerVolumetricDynamicScheme< 2, SparseSpaceType, LocalSpaceType > InnerVolumetricDynamicSchemeType2D;
    //typedef InnerVolumetricDynamicScheme< 3, SparseSpaceType, LocalSpaceType > InnerVolumetricDynamicSchemeType3D;


    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
    //typedef ResidualCriteria < SparseSpaceType, LocalSpaceType >::Pointer TResidual;
    //typedef DisplacementCriteria < SparseSpaceType, LocalSpaceType>::Pointer TDisplacement;


// 					;
    //********************************************************************
    //********************************************************************
    py::class_< ResidualBasedPredictorCorrectorBossakSchemeType, BaseSchemeType >
    (m, "ResidualBasedPredictorCorrectorBossakScheme")
    .def(py::init< double >());



    py::enum_<Constraint_Enforcement>(m, "Constraint_Enforcement")
    .value("Penalty_Methods", Penalty_Methods)
    .value("Lagrange_Multiplier_Methods", Lagrange_Multiplier_Methods)
    ;

    py::class_< ResidualBasedCentralDiferencesStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >, BaseSolvingStrategyType >
    (m, "ResidualBasedCentralDiferencesStrategy")
    .def(py::init< ModelPart&, Constraint_Enforcement, int, double, double, double, double,  bool, bool, bool, LinearSolverType::Pointer, BaseSchemeType::Pointer, BuilderAndSolverType::Pointer>())
    .def("Initialize", &ResidualBasedCentralDiferencesStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Initialize)
    .def("ComputeCriticalTime",  &ResidualBasedCentralDiferencesStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >::ComputeCriticalTime)
    .def("SetFractionDeltaTime", &ResidualBasedCentralDiferencesStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >::ChangeFractionDeltaTime)
    .def("SetConditionsFlag", &ResidualBasedCentralDiferencesStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >::ChangeContactConditions)
    .def("CalculateBoundaryContours", &ResidualBasedCentralDiferencesStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType >::RecalculateBoundaryContours)
    ;


    py::class_< ResidualBasedPredictorCorrectorBossakRotationSchemeType, BaseSchemeType >
    (m, "ResidualBasedPredictorCorrectorBossakRotationScheme")
    .def(py::init< double >());

    typedef ResidualBasedPredictorCorrectorRelaxationScheme< SparseSpaceType,
            LocalSpaceType > ResidualBasedPredictorCorrectorRelaxationSchemeType;

    py::class_< ResidualBasedPredictorCorrectorRelaxationSchemeType, BaseSchemeType >
    (m, "ResidualBasedPredictorCorrectorRelaxationScheme")
    .def(py::init< double, double >());

    py::class_< ResidualBasedNewmarkSchemeType, BaseSchemeType >
    (m, "ResidualBasedNewmarkScheme")
    .def(py::init< double >());

// 			py::class_< TestingSchemeType,
// 			bases< BaseSchemeType >,  boost::noncopyable >
// 					(
// 					"TestingScheme", py::init< >()
// 					);

    py::class_< MultiPhaseFlowCriteriaType, ConvergenceCriteriaBaseType >
    (m, "MultiPhaseFlowCriteria")
    .def(py::init<double, double >() );

//             py::class_< ResidualBasedMultiPhaseCriteriaType,
//             bases< ConvergenceCriteriaBaseType >, boost::noncopyable >
//             ("ResidualBasedMultiPhaseCriteria", py::init<double, double >() )
//             ;

            //py::class_ < MultiPhaseBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
            //( "MultiPhaseBuilderAndSolver", py::init<LinearSolverType::Pointer>() )
            //;



//            py::class_< ModalAnalysisBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
//                    (
//                     "ModalAnalysisBuilderAndSolver", py::init<LinearSolverType::Pointer>()
//                    );



//	py::class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeType,
//			bases< BaseSchemeType >,  boost::noncopyable >
//					(
//					"ResidualBasedPredictorCorrectorVelocityBossakScheme", py::init< double >()
//					);

    py::class_< VolumetricSchemeType2D, BaseSchemeType >
    (m, "VolumetricScheme2D")
    .def("CalculateCauchyStress",&VolumetricSchemeType2D::CalculateCauchyStress)
    ;

    py::class_< InnerVolumetricSchemeType2D, BaseSchemeType >
    (m, "InnerVolumetricScheme2D");

    py::class_< InnerVolumetricDynamicSchemeType2D, BaseSchemeType >
    (m, "InnerVolumetricDynamicScheme2D");

    py::class_< CompositSchemeType, BaseSchemeType >
    (m, "CompositScheme")
    .def(py::init< BaseSchemeType&, BaseSchemeType& >());

    py::class_< ResidualBasedArcLengthStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, BaseSolvingStrategyType >
    (m, "ResidualBasedArcLenghtStrategy")
    .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, unsigned int, unsigned int,double,bool, bool, bool,bool>())
    ;



    py::class_< ResidualBasedNewtonRaphsonLineSearchesStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, BaseSolvingStrategyType >
    (m, "ResidualBasedNewtonRaphsonLineSearchesStrategy")
    .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, unsigned int, unsigned int, double, double, double, double, bool, bool, bool, bool>() )
    ;




// 			py::class_<Residual_Displacement_Criteria<SparseSpaceType, LocalSpaceType >,
// 			         bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
// 			         boost::noncopyable >
// 			        ("ResidualDisplacementCriteria", py::init< double, double>() );
//
// 			py::class_<ResDisCriteria<SparseSpaceType, LocalSpaceType >,
// 			         bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
// 			         boost::noncopyable >
// 			        ("ResDisCriteria", py::init< TResidual,TDisplacement >());

    /*
    py::class_< ModalAnalysisBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable >
            (
             "ModalAnalysisBuilderAndSolver", py::init<LinearSolverType::Pointer>()
            );
    */
}
}  // namespace Python.

} // Namespace Kratos

