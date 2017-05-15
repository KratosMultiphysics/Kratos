// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
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
#include "custom_strategies/custom_strategies/line_search_contact_strategy.h"
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mesh_tying_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictionless_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictional_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_residual_contact_criteria.h"

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
    typedef LineSearchContactStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  LineSearchContactStrategyType;
    
    // Custom scheme types

    // Custom convergence criterion types
    typedef MeshTyingMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MeshTyingMortarConvergenceCriteriaType;
    typedef ALMFrictionlessMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionlessMortarConvergenceCriteriaType;
    typedef ALMFrictionalMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionalMortarConvergenceCriteriaType;
    typedef DisplacementLagrangeMultiplierContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierContactCriteriaType;
    typedef DisplacementLagrangeMultiplierResidualContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierResidualContactCriteriaType;
    
    // Custom builder and solvers types
    
    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************
            
    // Residual Based Newton Raphson Contact Strategy      
    class_< ResidualBasedNewtonRaphsonContactStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("ResidualBasedNewtonRaphsonContactStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
            .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
            .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
            .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
            .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
            ;
            
    // Residual Based Newton Raphson Contact Strategy      
    class_< LineSearchContactStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("LineSearchContactStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
            .def("SetMaxIterationNumber", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
            .def("GetMaxIterationNumber", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
            .def("SetKeepSystemConstantDuringIterations", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
            .def("GetKeepSystemConstantDuringIterations", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
            ;
            
    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************
            
    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Weighted residual values update
    class_< MeshTyingMortarConvergenceCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "MeshTyingMortarConvergenceCriteria", 
            init< >())
            .def(init< >())
            .def("SetEchoLevel", &MeshTyingMortarConvergenceCriteriaType::SetEchoLevel)
            ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless case)
    class_< ALMFrictionlessMortarConvergenceCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "ALMFrictionlessMortarConvergenceCriteria", 
            init< >())
            .def(init< >())
            .def("SetEchoLevel", &ALMFrictionlessMortarConvergenceCriteriaType::SetEchoLevel)
            ;
            
    // Dual set strategy for SSNM Convergence Criterion (frictional case)
    class_< ALMFrictionalMortarConvergenceCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "ALMFrictionalMortarConvergenceCriteria", 
            init< >())
            .def(init< >())
            .def("SetEchoLevel", &ALMFrictionalMortarConvergenceCriteriaType::SetEchoLevel)
            ;
            
    // Displacement and lagrange multiplier Convergence Criterion
    class_< DisplacementLagrangeMultiplierContactCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "DisplacementLagrangeMultiplierContactCriteria", 
            init< double, double, double, double >())
            .def("SetEchoLevel", &DisplacementLagrangeMultiplierContactCriteriaType::SetEchoLevel)
            ;
            
    // Displacement and lagrange multiplier residual Convergence Criterion
    class_< DisplacementLagrangeMultiplierResidualContactCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "DisplacementLagrangeMultiplierResidualContactCriteria", 
            init< double, double, double, double >())
            .def("SetEchoLevel", &DisplacementLagrangeMultiplierResidualContactCriteriaType::SetEchoLevel)
            ;
            
    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

}

}  // namespace Python.

} // Namespace Kratos

