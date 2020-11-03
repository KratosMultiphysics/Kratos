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
#include "custom_python/process_factory_utility.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/custom_strategies/line_search_contact_strategy.h"
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h"
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_mpc_contact_strategy.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mortar_and_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mesh_tying_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictionless_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/penalty_frictionless_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictionless_components_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/alm_frictional_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/penalty_frictional_mortar_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_frictional_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_mixed_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_mixed_frictional_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_residual_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_residual_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_lagrangemultiplier_residual_frictional_contact_criteria.h"
#include "custom_strategies/custom_convergencecriterias/contact_error_mesh_criteria.h"
#include "custom_strategies/custom_convergencecriterias/mpc_contact_criteria.h"

// Builders and solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_elimination_builder_and_solver_with_constraints.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

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
    typedef ResidualBasedNewtonRaphsonMPCContactStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedNewtonRaphsonMPCContactStrategyType;

    // Custom scheme types

    // Custom convergence criterion types
    typedef MortarAndConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MortarAndConvergenceCriteriaType;
    typedef MeshTyingMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > MeshTyingMortarConvergenceCriteriaType;
    typedef ALMFrictionlessMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionlessMortarConvergenceCriteriaType;
    typedef PenaltyFrictionlessMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > PenaltyFrictionlessMortarConvergenceCriteriaType;
    typedef ALMFrictionlessComponentsMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionlessComponentsMortarConvergenceCriteriaType;
    typedef ALMFrictionalMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ALMFrictionalMortarConvergenceCriteriaType;
    typedef PenaltyFrictionalMortarConvergenceCriteria< SparseSpaceType,  LocalSpaceType > PenaltyFrictionalMortarConvergenceCriteriaType;
    typedef DisplacementContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementContactCriteriaType;
    typedef DisplacementLagrangeMultiplierContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierContactCriteriaType;
    typedef DisplacementLagrangeMultiplierFrictionalContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierFrictionalContactCriteriaType;
    typedef DisplacementLagrangeMultiplierMixedContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierMixedContactCriteriaType;
    typedef DisplacementLagrangeMultiplierMixedFrictionalContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType;
    typedef DisplacementResidualContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementResidualContactCriteriaType;
    typedef DisplacementLagrangeMultiplierResidualContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierResidualContactCriteriaType;
    typedef DisplacementLagrangeMultiplierResidualFrictionalContactCriteria< SparseSpaceType,  LocalSpaceType > DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType;
    typedef ContactErrorMeshCriteria< SparseSpaceType,  LocalSpaceType > ContactErrorMeshCriteriaType;
    typedef MPCContactCriteria< SparseSpaceType,  LocalSpaceType > MPCContactCriteriaType;

    // Linear solvers

    // Custom builder and solvers types
    typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
    typedef ContactResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType, ResidualBasedBlockBuilderAndSolverType > ContactResidualBasedBlockBuilderAndSolverType;
    typedef ContactResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ContactResidualBasedEliminationBuilderAndSolverType;
    typedef ContactResidualBasedEliminationBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ContactResidualBasedEliminationBuilderAndSolverWithConstraintsType;

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Newton Raphson Contact Strategy
    py::class_< ResidualBasedNewtonRaphsonContactStrategyType,
        typename ResidualBasedNewtonRaphsonContactStrategyType::Pointer,
        BaseSolvingStrategyType  >  (m, "ResidualBasedNewtonRaphsonContactStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType>())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType, ProcessesListType>())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType>())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType, ProcessesListType>())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType>())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, ProcessesListType, ProcessesListType>())
        .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
        ;

    // Line search Contact Strategy
    py::class_< LineSearchContactStrategyType,
        typename LineSearchContactStrategyType::Pointer,
        BaseSolvingStrategyType  >(m, "LineSearchContactStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def("SetMaxIterationNumber", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
        ;

    // Residual Based Newton Raphson MPC Contact Strategy
    py::class_< ResidualBasedNewtonRaphsonMPCContactStrategyType, typename ResidualBasedNewtonRaphsonMPCContactStrategyType::Pointer, BaseSolvingStrategyType  >  (m, "ResidualBasedNewtonRaphsonMPCContactStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonMPCContactStrategyType::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonMPCContactStrategyType::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonMPCContactStrategyType::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonMPCContactStrategyType::GetKeepSystemConstantDuringIterations)
        ;

    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Custom mortar and criteria
    py::class_< MortarAndConvergenceCriteriaType, typename MortarAndConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "MortarAndConvergenceCriteria")
        .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer>())
        .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool>())
        .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool, ConditionNumberUtilityPointerType>())
        ;

    // Weighted residual values update
    py::class_< MeshTyingMortarConvergenceCriteriaType, typename MeshTyingMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "MeshTyingMortarConvergenceCriteria")
        .def(py::init< >())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless case)
    py::class_< ALMFrictionlessMortarConvergenceCriteriaType, typename ALMFrictionlessMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "ALMFrictionlessMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless penalty case)
    py::class_< PenaltyFrictionlessMortarConvergenceCriteriaType, typename PenaltyFrictionlessMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "PenaltyFrictionlessMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless components case)
    py::class_< ALMFrictionlessComponentsMortarConvergenceCriteriaType, typename ALMFrictionlessComponentsMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "ALMFrictionlessComponentsMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictional case)
    py::class_< ALMFrictionalMortarConvergenceCriteriaType, typename ALMFrictionalMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "ALMFrictionalMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        .def(py::init<bool, bool, bool, bool>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictional penalty case)
    py::class_< PenaltyFrictionalMortarConvergenceCriteriaType, typename PenaltyFrictionalMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "PenaltyFrictionalMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        .def(py::init<bool, bool, bool, bool>())
        ;

    // Displacement and lagrange multiplier Convergence Criterion
    py::class_< DisplacementContactCriteriaType, typename DisplacementContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double >())
        .def(py::init< double, double, double, double, bool >())
        ;

    // Displacement and lagrange multiplier Convergence Criterion
    py::class_< DisplacementLagrangeMultiplierContactCriteriaType, typename DisplacementLagrangeMultiplierContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementLagrangeMultiplierContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double, double, double >())
        .def(py::init< double, double, double, double, double, double, bool >())
        .def(py::init< double, double, double, double, double, double, bool, bool >())
        ;

    // Displacement and lagrange multiplier Convergence Criterion (frictional)
    py::class_< DisplacementLagrangeMultiplierFrictionalContactCriteriaType, typename DisplacementLagrangeMultiplierFrictionalContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementLagrangeMultiplierFrictionalContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double, bool >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double, bool, bool >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double, bool, bool, bool >())
        ;

    // Displacement and lagrange multiplier mixed Convergence Criterion
    py::class_< DisplacementLagrangeMultiplierMixedContactCriteriaType, typename DisplacementLagrangeMultiplierMixedContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementLagrangeMultiplierMixedContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double, double, double >())
        .def(py::init< double, double, double, double, double, double, bool >())
        .def(py::init< double, double, double, double, double, double, bool, bool >())
        ;

    // Displacement and lagrange multiplier mixed Convergence Criterion (frictional)
    py::class_< DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType, typename DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementLagrangeMultiplierMixedFrictionalContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double, bool >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double, bool, bool >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double, bool, bool, bool >())
        ;

    // Displacement residual Convergence Criterion
    py::class_< DisplacementResidualContactCriteriaType, typename DisplacementResidualContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementResidualContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double >())
        .def(py::init< double, double, double, double, bool >())
        ;

    // Displacement and lagrange multiplier residual Convergence Criterion
    py::class_< DisplacementLagrangeMultiplierResidualContactCriteriaType, typename DisplacementLagrangeMultiplierResidualContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementLagrangeMultiplierResidualContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double, double, double >())
        .def(py::init< double, double, double, double, double, double, bool >())
        .def(py::init< double, double, double, double, double, double, bool, bool >())
        ;

    // Displacement and lagrange multiplier residual Convergence Criterion (frictional)
    py::class_< DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType, typename DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "DisplacementLagrangeMultiplierResidualFrictionalContactCriteria")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double , bool >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double , bool, bool >())
        .def(py::init< double, double, double, double, double, double, double, double, double, double, double , bool, bool, bool >())
        ;

    // Error mesh Convergence Criterion
    py::class_< ContactErrorMeshCriteriaType, typename ContactErrorMeshCriteriaType::Pointer, ConvergenceCriteriaType >(m, "ContactErrorMeshCriteria")
        .def(py::init<Parameters>())
        ;

    // Contact convergence criteria
    py::class_< MPCContactCriteriaType, typename MPCContactCriteriaType::Pointer, ConvergenceCriteriaType > (m, "MPCContactCriteria")
        .def(py::init< >())
        ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    // Contact block builder and solver
    py::class_< ContactResidualBasedBlockBuilderAndSolverType, ContactResidualBasedBlockBuilderAndSolverType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedBlockBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer > ());

    // Contact elimination builder and solver
    py::class_< ContactResidualBasedEliminationBuilderAndSolverType, ContactResidualBasedEliminationBuilderAndSolverType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedEliminationBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer > ());

    // Contact elimination builder and sokver with constraints
    py::class_< ContactResidualBasedEliminationBuilderAndSolverWithConstraintsType, ContactResidualBasedEliminationBuilderAndSolverWithConstraintsType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedEliminationBuilderAndSolverWithConstraints")
    .def(py::init< LinearSolverType::Pointer > ());
}

}  // namespace Python.

} // Namespace Kratos

