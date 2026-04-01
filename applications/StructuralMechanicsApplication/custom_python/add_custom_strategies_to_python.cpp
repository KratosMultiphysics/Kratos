// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// Strategies
#include "custom_strategies/custom_strategies/eigensolver_strategy.hpp"
#include "custom_strategies/custom_strategies/modal_harmonic_analysis_strategy.hpp"
#include "custom_strategies/custom_strategies/direct_harmonic_analysis_strategy.hpp"
#include "custom_strategies/custom_strategies/formfinding_strategy.hpp"
#include "custom_strategies/custom_strategies/mechanical_explicit_strategy.hpp"
#include "custom_strategies/custom_strategies/prebuckling_strategy.hpp"

// Schemes
#include "custom_strategies/custom_schemes/residual_based_relaxation_scheme.hpp"
#include "custom_strategies/custom_schemes/explicit_central_differences_scheme.hpp"
#include "custom_strategies/custom_schemes/explicit_multi_stage_kim_scheme.hpp"
#include "custom_strategies/custom_schemes/eigensolver_dynamic_scheme.hpp"
#include "custom_strategies/custom_schemes/structural_mechanics_bossak_scheme.h"
#include "custom_strategies/custom_schemes/structural_mechanics_static_scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/residual_displacement_and_other_dof_criteria.h"
#include "custom_strategies/custom_convergencecriterias/error_mesh_criteria.h"

// Builders and solvers

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos::Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    using BaseBossakSchemeType = ResidualBasedBossakDisplacementScheme<SparseSpaceType, LocalSpaceType>;
    using BaseStaticSchemeType = ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ImplicitSolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
//     typedef BaseSolvingStrategyType::Pointer BaseSolvingStrategyPointer;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;
    typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

    // Custom strategy types
    typedef EigensolverStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > EigensolverStrategyType;
    typedef PrebucklingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > PrebucklingStrategyType;
    typedef ModalHarmonicAnalysisStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ModalHarmonicAnalysisStrategyType;
    typedef DirectHarmonicAnalysisStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > DirectHarmonicAnalysisStrategyType;
    typedef FormfindingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > FormfindingStrategyType;
    typedef MechanicalExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > MechanicalExplicitStrategyType;


    // Custom scheme types
    typedef ResidualBasedRelaxationScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRelaxationSchemeType;
    typedef EigensolverDynamicScheme< SparseSpaceType, LocalSpaceType > EigensolverDynamicSchemeType;
    typedef ExplicitCentralDifferencesScheme< SparseSpaceType, LocalSpaceType >  ExplicitCentralDifferencesSchemeType;
    typedef ExplicitMultiStageKimScheme< SparseSpaceType, LocalSpaceType >  ExplicitMultiStageKimSchemeType;
    using StructuralMechanicsBossakSchemeType = StructuralMechanicsBossakScheme<SparseSpaceType, LocalSpaceType>;
    using StructuralMechanicsStaticSchemeType = StructuralMechanicsStaticScheme<SparseSpaceType, LocalSpaceType>;

    // Custom convergence criterion types
    typedef ResidualDisplacementAndOtherDoFCriteria< SparseSpaceType,  LocalSpaceType > ResidualDisplacementAndOtherDoFCriteriaType;
    typedef ErrorMeshCriteria< SparseSpaceType,  LocalSpaceType > ErrorMeshCriteriaType;

    typedef std::complex<double> ComplexType;

    typedef UblasSpace<ComplexType,
                    boost::numeric::ublas::compressed_matrix<ComplexType>,
                    boost::numeric::ublas::vector<ComplexType>> ComplexSparseSpaceType;

    typedef UblasSpace<ComplexType,
                    boost::numeric::ublas::matrix<ComplexType>,
                    boost::numeric::ublas::vector<ComplexType>> ComplexLocalSpaceType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexLocalSpaceType> ComplexLinearSolverType;
    typedef ComplexLinearSolverType::Pointer ComplexLinearSolverPointer;

    // Custom builder and solvers types

    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Eigensolver Strategy
    py::class_< EigensolverStrategyType, typename EigensolverStrategyType::Pointer,BaseSolvingStrategyType >(m,"EigensolverStrategy")
        .def(py::init<ModelPart&,
             BaseSchemeType::Pointer,
             BuilderAndSolverPointer,
             double,
             double,
             bool,
             bool>(),
             py::arg("model_part"),
             py::arg("scheme"),
             py::arg("builder_and_solver"),
             py::arg("mass_matrix_diagonal_value"),
             py::arg("stiffness_matrix_diagonal_value"),
             py::arg("compute_modal_decomposition")=false,
             py::arg("normalize_eigenvectors_with_mass_matrix")=false)
            ;

    // Prebuckling Strategy
    py::class_< PrebucklingStrategyType, typename PrebucklingStrategyType::Pointer,BaseSolvingStrategyType >(m,"PrebucklingStrategy")
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverPointer, BuilderAndSolverPointer, ConvergenceCriteriaPointer, int, Parameters >())
        .def("GetSolutionFoundFlag", &PrebucklingStrategyType::GetSolutionFoundFlag)
        ;
    py::class_< FormfindingStrategyType,typename FormfindingStrategyType::Pointer, ResidualBasedNewtonRaphsonStrategyType >(m,"FormfindingStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaPointer, BuilderAndSolverPointer, ModelPart&, bool, const std::string&, Parameters, int, bool, bool, bool>())
        .def_static("WriteFormFoundMdpa", &FormfindingStrategyType::WriteFormFoundMdpa)
        ;


    py::class_< MechanicalExplicitStrategyType, typename MechanicalExplicitStrategyType::Pointer, BaseSolvingStrategyType >(m,"MechanicalExplicitStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, bool, bool, bool >())
        .def("SetInitializePerformedFlag", &MechanicalExplicitStrategyType::SetInitializePerformedFlag)
        .def("GetInitializePerformedFlag", &MechanicalExplicitStrategyType::GetInitializePerformedFlag)
        ;

    // harmonic Analysis Strategy
    py::class_< ModalHarmonicAnalysisStrategyType,typename ModalHarmonicAnalysisStrategyType::Pointer, BaseSolvingStrategyType >(m,"ModalHarmonicAnalysisStrategy")
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, BuilderAndSolverPointer, bool>() )
        .def("SetUseMaterialDampingFlag", &ModalHarmonicAnalysisStrategyType::SetUseMaterialDampingFlag)
        .def("GetUseMaterialDampingFlag", &ModalHarmonicAnalysisStrategyType::GetUseMaterialDampingFlag)
        ;

    py::class_<DirectHarmonicAnalysisStrategyType, typename DirectHarmonicAnalysisStrategyType::Pointer, BaseSolvingStrategyType>(m, "DirectHarmonicAnalysisStrategy")
        .def(py::init<ModelPart&,
                    BaseSchemeType::Pointer,
                    BuilderAndSolverPointer,
                    ComplexLinearSolverPointer,
                    double,
                    double,
                    double,
                    bool,
                    bool,
                    std::string,
                    std::string>(),
            py::arg("model_part"),
            py::arg("scheme"),
            py::arg("builder_and_solver"),
            py::arg("complex_linear_solver"),
            py::arg("mass_matrix_diagonal_value") = 1.0,
            py::arg("stiffness_matrix_diagonal_value") = 1.0,
            py::arg("damping_matrix_diagonal_value") = 1.0,
            py::arg("reform_dof_set_at_each_step") = false,
            py::arg("assemble_damping_matrix") = false,
            py::arg("real_load_sub_model_part") = "",
            py::arg("imaginary_load_sub_model_part") = "")
        .def("SetAssembleDampingMatrixFlag", &DirectHarmonicAnalysisStrategyType::SetAssembleDampingMatrixFlag)
        .def("GetAssembleDampingMatrixFlag", &DirectHarmonicAnalysisStrategyType::GetAssembleDampingMatrixFlag);

    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    // Residual Based Relaxation Scheme Type
    py::class_< ResidualBasedRelaxationSchemeType,typename ResidualBasedRelaxationSchemeType::Pointer, BaseSchemeType >(m,"ResidualBasedRelaxationScheme")
        .def(py::init< double , double >() )
        .def("Initialize", &ResidualBasedRelaxationScheme<SparseSpaceType, LocalSpaceType>::Initialize)
        ;

    // Eigensolver Scheme Type
    py::class_< EigensolverDynamicSchemeType,typename EigensolverDynamicSchemeType::Pointer, BaseSchemeType>(m,"EigensolverDynamicScheme")
        .def(py::init<>() )
        ;

    // Explicit Central Differences Scheme Type
    py::class_< ExplicitCentralDifferencesSchemeType,typename ExplicitCentralDifferencesSchemeType::Pointer, BaseSchemeType >(m,"ExplicitCentralDifferencesScheme")
        .def(py::init< const double, const double, const double>())
        .def(py::init< Parameters>())
        ;

    // Explicit Multi Stage Scheme Type
    py::class_< ExplicitMultiStageKimSchemeType,typename ExplicitMultiStageKimSchemeType::Pointer, BaseSchemeType >(m,"ExplicitMultiStageKimScheme")
        .def(py::init< const double>())
        .def(py::init< Parameters>())
        ;

    // Implicit Bossak Scheme Type
    py::class_<StructuralMechanicsBossakSchemeType, typename StructuralMechanicsBossakSchemeType::Pointer, BaseBossakSchemeType>(m, "StructuralMechanicsBossakScheme")
        .def(py::init<Parameters>())
        ;

    // Implicit Static Scheme Type
    py::class_<StructuralMechanicsStaticSchemeType, typename StructuralMechanicsStaticSchemeType::Pointer, BaseStaticSchemeType>(m, "StructuralMechanicsStaticScheme")
        .def(py::init<Parameters>())
        ;

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Displacement and other DoF residual Convergence Criterion
    py::class_< ResidualDisplacementAndOtherDoFCriteriaType,typename ResidualDisplacementAndOtherDoFCriteriaType::Pointer, ConvergenceCriteriaType >(m,"ResidualDisplacementAndOtherDoFCriteria")
        .def(py::init< double, double, const std::string >())
        .def(py::init< double, double>())
        ;

    // Error mesh Convergence Criterion
    py::class_< ErrorMeshCriteriaType, typename ErrorMeshCriteriaType::Pointer, ConvergenceCriteriaType >(m, "ErrorMeshCriteria")
        .def(py::init<Parameters>())
        ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************
}

}  // namespace Kratos::Python.

