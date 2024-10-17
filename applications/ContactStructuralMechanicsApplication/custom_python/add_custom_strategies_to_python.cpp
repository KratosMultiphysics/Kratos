// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
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

/* Strategies */
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_strategies/custom_strategies/line_search_contact_strategy.h"
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h"
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_mpc_contact_strategy.h"

/* Schemes */
#include "solving_strategies/schemes/scheme.h"

/* Convergence criterias */
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

/* Builders and solvers */
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solvers/contact_residualbased_elimination_builder_and_solver_with_constraints.h"

/* Linear solvers */
#include "linear_solvers/linear_solver.h"

namespace Kratos::Python
{
namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    using ProcessesListType = ProcessFactoryUtility::Pointer;
    using ConditionNumberUtilityPointerType = ConditionNumberUtility::Pointer;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;

    // Base types
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using BaseSolvingStrategyType = ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using ConvergenceCriteriaPointer = typename ConvergenceCriteriaType::Pointer;
    using BuilderAndSolverType = BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    // Custom strategy types
    using ResidualBasedNewtonRaphsonContactStrategyType = ResidualBasedNewtonRaphsonContactStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using LineSearchContactStrategyType = LineSearchContactStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using ResidualBasedNewtonRaphsonMPCContactStrategyType = ResidualBasedNewtonRaphsonMPCContactStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    // Custom scheme types

    // Custom convergence criterion types
    using MortarAndConvergenceCriteriaType = MortarAndConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using MeshTyingMortarConvergenceCriteriaType = MeshTyingMortarConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using ALMFrictionlessMortarConvergenceCriteriaType = ALMFrictionlessMortarConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using PenaltyFrictionlessMortarConvergenceCriteriaType = PenaltyFrictionlessMortarConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using ALMFrictionlessComponentsMortarConvergenceCriteriaType = ALMFrictionlessComponentsMortarConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using ALMFrictionalMortarConvergenceCriteriaType = ALMFrictionalMortarConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using PenaltyFrictionalMortarConvergenceCriteriaType = PenaltyFrictionalMortarConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementContactCriteriaType = DisplacementContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementLagrangeMultiplierContactCriteriaType = DisplacementLagrangeMultiplierContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementLagrangeMultiplierFrictionalContactCriteriaType = DisplacementLagrangeMultiplierFrictionalContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementLagrangeMultiplierMixedContactCriteriaType = DisplacementLagrangeMultiplierMixedContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementLagrangeMultiplierMixedFrictionalContactCriteriaType = DisplacementLagrangeMultiplierMixedFrictionalContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementResidualContactCriteriaType = DisplacementResidualContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementLagrangeMultiplierResidualContactCriteriaType = DisplacementLagrangeMultiplierResidualContactCriteria<SparseSpaceType, LocalSpaceType>;
    using DisplacementLagrangeMultiplierResidualFrictionalContactCriteriaType = DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<SparseSpaceType, LocalSpaceType>;
    using ContactErrorMeshCriteriaType = ContactErrorMeshCriteria<SparseSpaceType, LocalSpaceType>;
    using MPCContactCriteriaType = MPCContactCriteria<SparseSpaceType, LocalSpaceType>;

    // Linear solvers

    // Custom builder and solvers types
    using ResidualBasedBlockBuilderAndSolverType = ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using ContactResidualBasedBlockBuilderAndSolverType = ContactResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType, ResidualBasedBlockBuilderAndSolverType>;
    using ContactResidualBasedEliminationBuilderAndSolverType = ContactResidualBasedEliminationBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using ContactResidualBasedEliminationBuilderAndSolverWithConstraintsType = ContactResidualBasedEliminationBuilderAndSolverWithConstraints<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Newton Raphson Contact Strategy
    py::class_< ResidualBasedNewtonRaphsonContactStrategyType,
        typename ResidualBasedNewtonRaphsonContactStrategyType::Pointer,
        BaseSolvingStrategyType  >  (m, "ResidualBasedNewtonRaphsonContactStrategy")
        .def(py::init < ModelPart&, Parameters >())
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
        .def(py::init < ModelPart&, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters >())
        .def("SetMaxIterationNumber", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &LineSearchContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
        ;

    // Residual Based Newton Raphson MPC Contact Strategy
    py::class_< ResidualBasedNewtonRaphsonMPCContactStrategyType, typename ResidualBasedNewtonRaphsonMPCContactStrategyType::Pointer, BaseSolvingStrategyType  >  (m, "ResidualBasedNewtonRaphsonMPCContactStrategy")
        .def(py::init < ModelPart&, Parameters >())
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
        .def(py::init<Parameters>())
        .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer>())
        .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool>())
        .def(py::init<ConvergenceCriteriaPointer, ConvergenceCriteriaPointer, bool, ConditionNumberUtilityPointerType>())
        ;

    // Weighted residual values update
    py::class_< MeshTyingMortarConvergenceCriteriaType, typename MeshTyingMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "MeshTyingMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<Parameters>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless case)
    py::class_< ALMFrictionlessMortarConvergenceCriteriaType, typename ALMFrictionlessMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "ALMFrictionlessMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<Parameters>())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless penalty case)
    py::class_< PenaltyFrictionlessMortarConvergenceCriteriaType, typename PenaltyFrictionlessMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "PenaltyFrictionlessMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<Parameters>())
        .def(py::init<bool>())
        .def(py::init<bool, bool>())
        .def(py::init<bool, bool, bool>())
        ;

    // Dual set strategy for SSNM Convergence Criterion (frictionless components case)
    py::class_< ALMFrictionlessComponentsMortarConvergenceCriteriaType, typename ALMFrictionlessComponentsMortarConvergenceCriteriaType::Pointer,
        ConvergenceCriteriaType >
        (m, "ALMFrictionlessComponentsMortarConvergenceCriteria")
        .def(py::init< >())
        .def(py::init<Parameters>())
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
        .def(py::init<Parameters>())
        ;

    // Contact convergence criteria
    py::class_< MPCContactCriteriaType, typename MPCContactCriteriaType::Pointer, ConvergenceCriteriaType > (m, "MPCContactCriteria")
        .def(py::init< >())
        .def(py::init<Parameters>())
        ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    // Contact block builder and solver
    py::class_< ContactResidualBasedBlockBuilderAndSolverType, ContactResidualBasedBlockBuilderAndSolverType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedBlockBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer > ())
    .def(py::init< LinearSolverType::Pointer, Parameters > ())
    ;

    // Contact elimination builder and solver
    py::class_< ContactResidualBasedEliminationBuilderAndSolverType, ContactResidualBasedEliminationBuilderAndSolverType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedEliminationBuilderAndSolver")
    .def(py::init< LinearSolverType::Pointer > ())
    .def(py::init< LinearSolverType::Pointer, Parameters > ())
    ;

    // Contact elimination builder and sokver with constraints
    py::class_< ContactResidualBasedEliminationBuilderAndSolverWithConstraintsType, ContactResidualBasedEliminationBuilderAndSolverWithConstraintsType::Pointer, BuilderAndSolverType > (m, "ContactResidualBasedEliminationBuilderAndSolverWithConstraints")
    .def(py::init< LinearSolverType::Pointer > ())
    .def(py::init< LinearSolverType::Pointer, Parameters > ())
    ;
}

} // Namespace Kratos::Python.

