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

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// Strategies
#include "custom_strategies/custom_strategies/residual_based_arc_length_strategy.hpp"
#include "custom_strategies/custom_strategies/eigensolver_strategy.hpp"
#include "custom_strategies/custom_strategies/harmonic_analysis_strategy.hpp"
#include "custom_strategies/custom_strategies/formfinding_updated_reference_strategy.hpp"
#include "custom_strategies/custom_strategies/mechanical_explicit_strategy.hpp"
#include "custom_strategies/custom_strategies/newton_raphson_with_hydrostatic_load_strategy.hpp"

// Schemes
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"
#include "custom_strategies/custom_schemes/explicit_central_differences_scheme.hpp"
#include "custom_strategies/custom_schemes/eigensolver_dynamic_scheme.hpp"
#include "custom_response_functions/adjoint_schemes/adjoint_structural_static_scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/displacement_and_other_dof_criteria.h"
#include "custom_strategies/custom_convergencecriterias/residual_displacement_and_other_dof_criteria.h"
#include "custom_strategies/custom_convergencecriterias/error_mesh_criteria.h"

// Builders and solvers

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{
using namespace pybind11;

void AddCustomStrategiesToPython(pybind11::module &m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef LinearSolverType::Pointer LinearSolverPointer;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> BaseSolvingStrategyType;
    //     typedef BaseSolvingStrategyType::Pointer BaseSolvingStrategyPointer;
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
    typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;
    typedef ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedNewtonRaphsonStrategyType;

    // Custom strategy types
    // typedef ResidualBasedArcLengthStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedArcLengthStrategyType;
    typedef EigensolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> EigensolverStrategyType;
    typedef HarmonicAnalysisStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> HarmonicAnalysisStrategyType;
    typedef FormfindingUpdatedReferenceStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> FormfindingUpdatedReferenceStrategyType;
    typedef MechanicalExplicitStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> MechanicalExplicitStrategyType;
    typedef NewtonRaphsonWithHydrostaticLoadStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> NewtonRaphsonWithHydrostaticLoadStrategyType;

    // Custom scheme types
    typedef ResidualBasedRelaxationScheme<SparseSpaceType, LocalSpaceType> ResidualBasedRelaxationSchemeType;
    typedef EigensolverDynamicScheme<SparseSpaceType, LocalSpaceType> EigensolverDynamicSchemeType;
    typedef ExplicitCentralDifferencesScheme<SparseSpaceType, LocalSpaceType> ExplicitCentralDifferencesSchemeType;
    typedef AdjointStructuralStaticScheme<SparseSpaceType, LocalSpaceType> AdjointStructuralStaticSchemeType;

    // Custom convergence criterion types
    typedef DisplacementAndOtherDoFCriteria<SparseSpaceType, LocalSpaceType> DisplacementAndOtherDoFCriteriaType;
    typedef ResidualDisplacementAndOtherDoFCriteria<SparseSpaceType, LocalSpaceType> ResidualDisplacementAndOtherDoFCriteriaType;
    typedef ErrorMeshCriteria<SparseSpaceType, LocalSpaceType> ErrorMeshCriteriaType;

    // Custom builder and solvers types

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Arc Length Strategy
    // Currently not woking
    // class_< ResidualBasedArcLengthStrategyType,typename ResidualBasedArcLengthStrategyType::Pointer, BaseSolvingStrategyType >(m,"ResidualBasedArcLengthStrategy")
    // .def(init<ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer,
    //                                                             unsigned int, unsigned int, unsigned int,long double,bool, bool, bool>() )
    //        ;

    // Eigensolver Strategy
    class_<EigensolverStrategyType, typename EigensolverStrategyType::Pointer, BaseSolvingStrategyType>(m, "EigensolverStrategy")
        .def(init<ModelPart &, BaseSchemeType::Pointer, BuilderAndSolverPointer>());

    class_<FormfindingUpdatedReferenceStrategyType, typename FormfindingUpdatedReferenceStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType>(m, "FormfindingUpdatedReferenceStrategy")
        .def(init<ModelPart &, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, int, bool, bool, bool, bool, bool>())
        .def(init<ModelPart &, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, BuilderAndSolverPointer, int, bool, bool, bool, bool, bool>())
        .def("SetMaxIterationNumber", &FormfindingUpdatedReferenceStrategyType::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &FormfindingUpdatedReferenceStrategyType::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations", &FormfindingUpdatedReferenceStrategyType::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations", &FormfindingUpdatedReferenceStrategyType::GetKeepSystemConstantDuringIterations)
        .def("SetInitializePerformedFlag", &FormfindingUpdatedReferenceStrategyType::SetInitializePerformedFlag)
        .def("GetInitializePerformedFlag", &FormfindingUpdatedReferenceStrategyType::GetInitializePerformedFlag);

    //hydrostatic Analysis Strategy

    class_<NewtonRaphsonWithHydrostaticLoadStrategyType, typename NewtonRaphsonWithHydrostaticLoadStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType>(m, "NewtonRaphsonWithHydrostaticLoadStrategy")
        .def(init<ModelPart &, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, int, bool, bool, bool, bool>())
        .def(init<ModelPart &, BaseSchemeType::Pointer, LinearSolverPointer, ConvergenceCriteriaPointer, BuilderAndSolverPointer, int, bool, bool, bool, bool>());

    class_<MechanicalExplicitStrategyType, typename MechanicalExplicitStrategyType::Pointer, BaseSolvingStrategyType>(m, "MechanicalExplicitStrategy")
        .def(init<ModelPart &, BaseSchemeType::Pointer, bool, bool, bool>())
        .def("SetInitializePerformedFlag", &MechanicalExplicitStrategyType::SetInitializePerformedFlag)
        .def("GetInitializePerformedFlag", &MechanicalExplicitStrategyType::GetInitializePerformedFlag);

    // harmonic Analysis Strategy
    class_<HarmonicAnalysisStrategyType, typename HarmonicAnalysisStrategyType::Pointer, BaseSolvingStrategyType>(m, "HarmonicAnalysisStrategy")
        .def(init<ModelPart &, BaseSchemeType::Pointer, BuilderAndSolverPointer, bool>())
        .def("SetUseMaterialDampingFlag", &HarmonicAnalysisStrategyType::SetUseMaterialDampingFlag)
        .def("GetUseMaterialDampingFlag", &HarmonicAnalysisStrategyType::GetUseMaterialDampingFlag);

    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    // Residual Based Relaxation Scheme Type
    class_<ResidualBasedRelaxationSchemeType, typename ResidualBasedRelaxationSchemeType::Pointer, BaseSchemeType>(m, "ResidualBasedRelaxationScheme")
        .def(init<double, double>())
        .def("Initialize", &ResidualBasedRelaxationScheme<SparseSpaceType, LocalSpaceType>::Initialize);

    // Eigensolver Scheme Type
    class_<EigensolverDynamicSchemeType, typename EigensolverDynamicSchemeType::Pointer, BaseSchemeType>(m, "EigensolverDynamicScheme")
        .def(init<>());

    // Explicit Central Differences Scheme Type
    class_<ExplicitCentralDifferencesSchemeType, typename ExplicitCentralDifferencesSchemeType::Pointer, BaseSchemeType>(m, "ExplicitCentralDifferencesScheme")
        .def(init<const double, const double, const double>())
        .def(init<Parameters>());

    class_<AdjointStructuralStaticSchemeType, AdjointStructuralStaticSchemeType::Pointer, BaseSchemeType>(m, "AdjointStructuralStaticScheme")
        .def(init<Parameters, AdjointStructuralResponseFunction::Pointer>());

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Displacement and other DoF Convergence Criterion
    class_<DisplacementAndOtherDoFCriteriaType, typename DisplacementAndOtherDoFCriteriaType::Pointer, ConvergenceCriteriaType>(m, "DisplacementAndOtherDoFCriteria")
        .def(init<double, double, std::string>())
        .def(init<double, double>());

    // Displacement and other DoF residual Convergence Criterion
    class_<ResidualDisplacementAndOtherDoFCriteriaType, typename ResidualDisplacementAndOtherDoFCriteriaType::Pointer, ConvergenceCriteriaType>(m, "ResidualDisplacementAndOtherDoFCriteria")
        .def(init<double, double, std::string>())
        .def(init<double, double>());

    // Error mesh Convergence Criterion
    class_<ErrorMeshCriteriaType, typename ErrorMeshCriteriaType::Pointer, ConvergenceCriteriaType>(m, "ErrorMeshCriteria")
        .def(init<Parameters>());

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************
}

} // namespace Python.

} // Namespace Kratos
