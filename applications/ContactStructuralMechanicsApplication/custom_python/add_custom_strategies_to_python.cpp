// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix
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
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/custom_schemes/residual_based_incremental_update_static_contact_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_incremental_update_static_ALM_contact_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_bossak_displacement_contact_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_bossak_displacement_ALM_contact_scheme.hpp"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mortar_criteria.h"

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
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
        
    // Custom strategy types
    typedef ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedNewtonRaphsonContactStrategyType;
    
    // Custom scheme types
    typedef ResidualBasedIncrementalUpdateStaticContactScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticContactSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticALMContactScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticALMContactSchemeType;
    typedef ResidualBasedBossakDisplacementContactScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementContactSchemeType;
    typedef ResidualBasedBossakDisplacementALMContactScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementALMContactSchemeType;

    // Custom convergence criterion types
    typedef MortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MortarConvergenceCriteriaType;
    
    // Custom builder and solvers types
    
    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
            
    // Residual Based Newton Raphson Contact Strategy      
    class_< ResidualBasedNewtonRaphsonContactStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("ResidualBasedNewtonRaphsonContactStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, double, unsigned int >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, double, unsigned int >())
            .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
            .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
            .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
            .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
            ;
            
    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    // Residual Based Incremental Update Static Contact Scheme Type
    class_< ResidualBasedIncrementalUpdateStaticContactSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
            "ResidualBasedIncrementalUpdateStaticContactScheme", init< >()
            );
            
    // Residual Based Incremental Update Static Contact Scheme Type
    class_< ResidualBasedIncrementalUpdateStaticALMContactSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
            "ResidualBasedIncrementalUpdateStaticALMContactScheme", init< >()
            );
            
    // Residual Based Bossak Scheme Type
    class_< ResidualBasedBossakDisplacementContactSchemeType,
    bases< BaseSchemeType >,  boost::noncopyable >
    (
        "ResidualBasedBossakDisplacementContactScheme", init< double >() )
        .def("Initialize", &ResidualBasedBossakDisplacementContactScheme<SparseSpaceType, LocalSpaceType>::Initialize)
    ;
    
    // Residual Based Bossak Scheme Type
    class_< ResidualBasedBossakDisplacementALMContactSchemeType,
    bases< BaseSchemeType >,  boost::noncopyable >
    (
        "ResidualBasedBossakDisplacementALMContactScheme", init< double >() )
        .def("Initialize", &ResidualBasedBossakDisplacementALMContactScheme<SparseSpaceType, LocalSpaceType>::Initialize)
    ;
     
    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Displacement Convergence Criterion
    class_< MortarConvergenceCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "MortarConvergenceCriteria", 
            init< >())
            .def(init< >())
            .def("SetEchoLevel", &MortarConvergenceCriteriaType::SetEchoLevel)
            ;
            
    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

}

}  // namespace Python.

} // Namespace Kratos

