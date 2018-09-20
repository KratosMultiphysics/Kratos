// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_utilities/process_factory_utility.h"
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
#include "custom_strategies/custom_convergencecriterias/mortar_and_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mesh_tying_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictionless_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictionless_components_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictional_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_frictional_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_mixed_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_mixed_frictional_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_residual_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_residual_frictional_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/contact_error_mesh_criteria.h"

// Builders and solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_block_builder_and_solver.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{
using namespace pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef ProcessFactoryUtility::Pointer ProcessesListType;
    typedef ConditionNumberUtility::Pointer ConditionNumberUtilityPointerType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

    // Custom strategy types
    typedef ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedNewtonRaphsonContactStrategyType;
    typedef LineSearchContactStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  LineSearchContactStrategyType;

    // Custom scheme types

    // Custom convergence criterion types
    typedef MortarAndConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MortarAndConvergenceCriteriaType;
    typedef MeshTyingMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MeshTyingMortarConvergenceCriteriaType;
    typedef ALMFrictionlessMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionlessMortarConvergenceCriteriaType;
    typedef ALMFrictionlessComponentsMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionlessComponentsMortarConvergenceCriteriaType;
    typedef ALMFrictionalMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionalMortarConvergenceCriteriaType;
    typedef DisplacementLagrangeMultiplierContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierContactCriteriaType;
    typedef DisplacementLagrangeMultiplierFrictionalContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierFrictionalContactCriteriaType;
    typedef DisplacementLagrangeMultiplierMixedContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierMixedContactCriteriaType;
    typedef DisplacementLagrangeMultiplierMixedFrictionalContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType;
    typedef DisplacementLagrangeMultiplierResidualContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierResidualContactCriteriaType;
    typedef DisplacementLagrangeMultiplierResidualFrictionalContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType;
    typedef ContactErrorMeshCriteria< SparseSpaceType,  LocalSpaceType > ContactErrorMeshCriteriaType;

    // Linear solvers

    // Custom builder and solvers types
    typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsType;
    typedef ContactResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, ResidualBasedBlockBuilderAndSolverType > ContactResidualBasedBlockBuilderAndSolverType;
    typedef ContactResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, ResidualBasedBlockBuilderAndSolverWithConstraintsType > ContactResidualBasedBlockBuilderAndSolverWithConstraintsType;

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Newton Raphson Contact Strategy
    class_< ResidualBasedNewtonRaphsonContactStrategyType,
            typename ResidualBasedNewtonRaphsonContactStrategyType::Pointer,
            BaseSolvingStrategyType  >  (m, "ResidualBasedNewtonRaphsonContactStrategy")
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType>())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType, ProcessesListType>())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType>())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType, ProcessesListType>())
            .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
            .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
            .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
            .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
            ;

    // Line search Contact Strategy
    class_< LineSearchContactStrategyType,
            typename LineSearchContactStrategyType::Pointer,
            BaseSolvingStrategyType  >(m, "LineSearchContactStrategy")
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
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

    // Custom mortar and criteria
    class_< MortarAndConvergenceCriteriaType, typename MortarAndConvergenceCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "MortarAndConvergenceCriteria")
            .def(init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer>())
            .def(init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool>())
            .def(init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool, ConditionNumberUtilityPointerType>())
            ;

    // Weighted residual values update
    class_< MeshTyingMortarConvergenceCriteriaType, typename MeshTyingMortarConvergenceCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "MeshTyingMortarConvergenceCriteria")
            .def(init< >())
            ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless case)
    class_< ALMFrictionlessMortarConvergenceCriteriaType, typename ALMFrictionlessMortarConvergenceCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "ALMFrictionlessMortarConvergenceCriteria")
            .def(init< >())
            .def(init<bool>())
            .def(init<bool, bool>())
            ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless components case)
    class_< ALMFrictionlessComponentsMortarConvergenceCriteriaType, typename ALMFrictionlessComponentsMortarConvergenceCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "ALMFrictionlessComponentsMortarConvergenceCriteria")
            .def(init< >())
            .def(init<bool>())
            .def(init<bool, bool>())
            ;

    // Dual set strategy for SSNM Convergence Criterion (frictional case)
    class_< ALMFrictionalMortarConvergenceCriteriaType, typename ALMFrictionalMortarConvergenceCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "ALMFrictionalMortarConvergenceCriteria")
            .def(init< >())
            .def(init<bool>())
            .def(init<bool, bool>())
            ;

    // Displacement and lagrange multiplier Convergence Criterion
    class_< DisplacementLagrangeMultiplierContactCriteriaType, typename DisplacementLagrangeMultiplierContactCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "DisplacementLagrangeMultiplierContactCriteria")
            .def(init<>())
            .def(init<Parameters>())
            .def(init< double, double, double, double >())
            .def(init< double, double, double, double, bool >())
            .def(init< double, double, double, double, bool, bool >())
            ;

    // Displacement and lagrange multiplier Convergence Criterion (frictional)
    class_< DisplacementLagrangeMultiplierFrictionalContactCriteriaType, typename DisplacementLagrangeMultiplierFrictionalContactCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "DisplacementLagrangeMultiplierFrictionalContactCriteria")
            .def(init<>())
            .def(init<Parameters>())
            .def(init< double, double, double, double, double, double  >())
            .def(init< double, double, double, double, double, double , bool >())
            .def(init< double, double, double, double, double, double , bool, bool >())
            ;
            
    // Displacement and lagrange multiplier mixed Convergence Criterion
    class_< DisplacementLagrangeMultiplierMixedContactCriteriaType, typename DisplacementLagrangeMultiplierMixedContactCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "DisplacementLagrangeMultiplierMixedContactCriteria")
            .def(init<>())
            .def(init<Parameters>())
            .def(init< double, double, double, double >())
            .def(init< double, double, double, double, bool >())
            .def(init< double, double, double, double, bool, bool >())
            ;
  
    // Displacement and lagrange multiplier mixed Convergence Criterion (frictional)
    class_< DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType, typename DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "DisplacementLagrangeMultiplierMixedFrictionalContactCriteria")
            .def(init<>())
            .def(init<Parameters>())
            .def(init< double, double, double, double, double, double >())
            .def(init< double, double, double, double, double, double, bool >())
            .def(init< double, double, double, double, double, double, bool, bool >())
            ;
            
    // Displacement and lagrange multiplier residual Convergence Criterion
    class_< DisplacementLagrangeMultiplierResidualContactCriteriaType, typename DisplacementLagrangeMultiplierResidualContactCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "DisplacementLagrangeMultiplierResidualContactCriteria")
            .def(init<>())
            .def(init<Parameters>())
            .def(init< double, double, double, double >())
            .def(init< double, double, double, double, bool >())
            .def(init< double, double, double, double, bool, bool >())
            ;

    // Displacement and lagrange multiplier residual Convergence Criterion (frictional)
    class_< DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType, typename DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType::Pointer,
            ConvergenceCriteriaType >
            (m, "DisplacementLagrangeMultiplierResidualFrictionalContactCriteria")
            .def(init<>())
            .def(init<Parameters>())
            .def(init< double, double, double, double, double, double >())
            .def(init< double, double, double, double, double, double , bool >())
            .def(init< double, double, double, double, double, double , bool, bool >())
            ;
            
    // Error mesh Convergence Criterion
    class_< ContactErrorMeshCriteriaType, typename ContactErrorMeshCriteriaType::Pointer, ConvergenceCriteriaType >(m, "ContactErrorMeshCriteria")
    .def(init<Parameters>())
    ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    // Contact block builder and solver
    class_< ContactResidualBasedBlockBuilderAndSolverType, ContactResidualBasedBlockBuilderAndSolverType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedBlockBuilderAndSolver")
    .def(init< LinearSolverType::Pointer > ());

    // Contact block buiklder and sokver with constraints
    class_< ContactResidualBasedBlockBuilderAndSolverWithConstraintsType, ContactResidualBasedBlockBuilderAndSolverWithConstraintsType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedBlockBuilderAndSolverWithConstraints")
    .def(init< LinearSolverType::Pointer > ());
}

}  // namespace Python.

} // Namespace Kratos

