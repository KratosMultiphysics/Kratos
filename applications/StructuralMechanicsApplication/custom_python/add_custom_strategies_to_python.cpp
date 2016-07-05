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
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_incremental_update_static_contact_scheme.hpp"
#include "custom_strategies/custom_strategies/residual_based_arc_length_strategy.hpp"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mortar_criteria.h"

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

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    
    // Custom strategy types
     typedef ResidualBasedArcLengthStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedArcLengthStrategyType;
     
    // Custom builder_and_solver types

    // Custom scheme types
    typedef ResidualBasedRelaxationScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRelaxationSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticContactScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticContactSchemeType;

    // Custom convergence criterion types
    typedef MortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MortarConvergenceCriteriaType;
    
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

    // Residual Based Incremental Update Static Contact Scheme Type
    class_< ResidualBasedIncrementalUpdateStaticContactSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
            "ResidualBasedIncrementalUpdateStaticContactScheme", init< >()
            );
     
    // Residual Based Arc Length Strategy      
    class_< ResidualBasedArcLengthStrategyType,
            bases< BaseSolvingStrategyType >,  boost::noncopyable >
            (
                "ResidualBasedArcLengthStrategy", init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                                                                unsigned int, unsigned int, unsigned int,long double,bool, bool, bool>() )
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
            
}

}  // namespace Python.

} // Namespace Kratos

