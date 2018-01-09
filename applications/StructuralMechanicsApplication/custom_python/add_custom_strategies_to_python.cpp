// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
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

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/custom_strategies/residual_based_arc_length_strategy.hpp"
#include "custom_strategies/custom_strategies/eigensolver_strategy.hpp"
#include "custom_strategies/custom_strategies/harmonic_analysis_strategy.hpp"
#include "custom_strategies/custom_strategies/formfinding_updated_reference_strategy.hpp"

// Schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"
#include "custom_strategies/custom_schemes/eigensolver_dynamic_scheme.hpp"

// Builder and solvers
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solver/residualbased_block_builder_and_solver_with_mpc.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_and_other_dof_criteria.h"
#include "custom_strategies/custom_convergencecriterias/residual_displacement_and_other_dof_criteria.h"

// Builders and solvers

// Linear solvers
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
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef LinearSolverType::Pointer LinearSolverPointer;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
//     typedef BaseSolvingStrategyType::Pointer BaseSolvingStrategyPointer;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;
    
    // Custom strategy types
    typedef ResidualBasedArcLengthStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedArcLengthStrategyType;
    typedef EigensolverStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > EigensolverStrategyType;
    typedef HarmonicAnalysisStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > HarmonicAnalysisStrategyType;
    typedef FormfindingUpdatedReferenceStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > FormfindingUpdatedReferenceStrategyType;

    // Custom scheme types
    typedef ResidualBasedRelaxationScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRelaxationSchemeType;
    typedef EigensolverDynamicScheme< SparseSpaceType, LocalSpaceType > EigensolverDynamicSchemeType;
    
    // Custom convergence criterion types
    typedef DisplacementAndOtherDoFCriteria< SparseSpaceType,  LocalSpaceType > DisplacementAndOtherDoFCriteriaType;
    typedef ResidualDisplacementAndOtherDoFCriteria< SparseSpaceType,  LocalSpaceType > ResidualDisplacementAndOtherDoFCriteriaType;

    // Custom builder and solvers types

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
    
    // Residual Based Arc Length Strategy      
    class_< ResidualBasedArcLengthStrategyType, bases< BaseSolvingStrategyType >,  boost::noncopyable >
            (
                "ResidualBasedArcLengthStrategy", init<ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer,
                                                                unsigned int, unsigned int, unsigned int,long double,bool, bool, bool>() )
            ;

    // Eigensolver Strategy
    class_< EigensolverStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
            (
                "EigensolverStrategy", init<ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverPointer>() )
            ;
             

    class_< FormfindingUpdatedReferenceStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
        ("FormfindingUpdatedReferenceStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, int, bool, bool, bool >())
        .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, BuilderAndSolverPointer, int, bool, bool, bool >())
        .def("SetMaxIterationNumber", &FormfindingUpdatedReferenceStrategyType::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &FormfindingUpdatedReferenceStrategyType::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &FormfindingUpdatedReferenceStrategyType::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &FormfindingUpdatedReferenceStrategyType::GetKeepSystemConstantDuringIterations)
        .def("SetInitializePerformedFlag", &FormfindingUpdatedReferenceStrategyType::SetInitializePerformedFlag)
        .def("GetInitializePerformedFlag", &FormfindingUpdatedReferenceStrategyType::GetInitializePerformedFlag)
        ;

    // harmonic Analysis Strategy
    class_< HarmonicAnalysisStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
            (
                "HarmonicAnalysisStrategy", init<ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverPointer, bool>() )
            // .def(init < ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverPointer >())
            .def("SetUseMaterialDampingFlag", &HarmonicAnalysisStrategyType::SetUseMaterialDampingFlag)
            .def("GetUseMaterialDampingFlag", &HarmonicAnalysisStrategyType::GetUseMaterialDampingFlag)
            ;

    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************
    
    // Residual Based Relaxation Scheme Type
    class_< ResidualBasedRelaxationSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedRelaxationScheme", init< double , double >() )
            .def("Initialize", &ResidualBasedRelaxationScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;

    // Eigensolver Scheme Type
    class_< EigensolverDynamicSchemeType,
            EigensolverDynamicSchemeType::Pointer, bases< BaseSchemeType >, boost::noncopyable >
            (
                "EigensolverDynamicScheme", init<>() )
            ;

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************
            
    // Displacement and other DoF Convergence Criterion
    class_< DisplacementAndOtherDoFCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "DisplacementAndOtherDoFCriteria", 
            init< double, double, std::string >())
            .def(init< double, double>())
            ;
            
    // Displacement and other DoF residual Convergence Criterion
    class_< ResidualDisplacementAndOtherDoFCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "ResidualDisplacementAndOtherDoFCriteria", 
            init< double, double, std::string >())
            .def(init< double, double>())
            ;
            
    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************
    class_< ResidualBasedBlockBuilderAndSolverWithMpc< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                bases< ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >,
                boost::noncopyable >
                ("ResidualBasedBlockBuilderAndSolverWithMpc", init<LinearSolverType::Pointer>());
}

}  // namespace Python.

} // Namespace Kratos

