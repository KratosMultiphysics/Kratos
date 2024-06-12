//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:    Hubert Balcerzak
//                   Riccardo Rossi
//

#if !defined(CSHARP_WRAPPER_APPLICATION_KRATOS_INTERNALS_H_INCLUDED )
#define  CSHARP_WRAPPER_APPLICATION_KRATOS_INTERNALS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kernel.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application.h"
#include "convection_diffusion_application.h"
#include "constitutive_laws_application.h"
#include "rom_application.h"
#include "includes/model_part_io.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "custom_strategies/strategies/residualbased_convdiff_strategy.h"
#include "custom_strategies/strategies/residualbased_convdiff_strategy_nonlinear.h"
#include "custom_strategies/strategies/residualbased_eulerian_convdiff_strategy.h"
#include "custom_strategies/rom_builder_and_solver.h"
#include "includes/constitutive_law.h"
#include "factories/linear_solver_factory.h"

namespace CSharpKratosWrapper {


    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrix;
    typedef boost::numeric::ublas::vector<double> Vector;
    typedef Kratos::UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef boost::numeric::ublas::matrix<double> Matrix;
    typedef Kratos::UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef Kratos::UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef Kratos::LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    typedef Kratos::ROMBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ROMBuilderAndSolverType;
    typedef Kratos::ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;

    typedef Kratos::ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;

    typedef Kratos::ResidualCriteria<SparseSpaceType, LocalSpaceType > ResidualCriteriaType;

    typedef Kratos::ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;
    typedef Kratos::ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;
    //ConvectionDiffusion Strategies
    typedef Kratos::ResidualBasedConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedConvectionDiffusionStrategyType;
    typedef Kratos::ResidualBasedConvectionDiffusionStrategyNonLinear< SparseSpaceType, LocalSpaceType, LinearSolverType >ResidualBasedConvectionDiffusionStrategyNonLinearType;
    typedef Kratos::ResidualBasedEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >ResidualBasedEulerianConvectionDiffusionStrategyType;


    /// The definition of the linear solver factory type
    typedef Kratos::LinearSolverFactory<SparseSpaceType,LocalSpaceType> LinearSolverFactoryType;

    class KratosInternals {

    public:
        void initInternals();
        void initInternalsRom();
        void initInternalsConvDiff();
        void initModelPart();
        void loadMDPA(const std::string& rMDPAFilePath);
        void loadSettingsParameters(const std::string& rJSONFilePath);
        void loadRomParameters(const std::string& rJSONFilePath);
        void loadHRomParameters(const std::string& rJSONFilePath);
        void initRomGeometry();
        void initHRomWeight();
        void initDofs();
        void initProperties();
        void initSolver();
        void initSolverRom();
        void initSolverConvDiff();
        void solve();
        Kratos::ModelPart& GetMainModelPart();
        Kratos::Parameters GetSettings();

    private:
        Kratos::Kernel mKernel;
        Kratos::KratosStructuralMechanicsApplication mApplication;
        Kratos::KratosConvectionDiffusionApplication mConvectionDiffusion;
        Kratos::KratosConstitutiveLawsApplication mConstitutiveLaws;
        Kratos::KratosRomApplication mRomApplication;
        std::string mModelpartName;
        Kratos::Model mModel;
        Kratos::Parameters mSettingsParameters;
        Kratos::Parameters mRomParameters;
        Kratos::Parameters mHRomParameters;
        ResidualBasedLinearStrategyType::Pointer pmLinearStrategy;
        ResidualBasedNewtonRaphsonStrategyType::Pointer pmStrategy;
        ResidualBasedConvectionDiffusionStrategyType::Pointer pmConvStrategy;
        ResidualBasedConvectionDiffusionStrategyNonLinearType::Pointer pmNonLinearConvStrategy;
        bool mRomEnable;
        bool mConvDiffEnabled;

        Kratos::Parameters GetDefaultParameters();
        Kratos::Parameters GetDefaultRomParameters();
        Kratos::Parameters GetDefaultHRomParameters();
    };
}

#endif	/* CSHARP_WRAPPER_APPLICATION_KRATOS_INTERNALS_H_INCLUDED */

