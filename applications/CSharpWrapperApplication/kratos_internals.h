#include "includes/kernel.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application.h"
#include "includes/model_part_io.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
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
    typedef Kratos::ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;

    typedef Kratos::ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;

    typedef Kratos::ResidualCriteria<SparseSpaceType, LocalSpaceType > ResidualCriteriaType;

    typedef Kratos::ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;

    /// The definition of the linear solver factory type
    typedef Kratos::LinearSolverFactory<SparseSpaceType,LocalSpaceType> LinearSolverFactoryType;

    class KratosInternals {

    public:
        void initInternals();
        void initModelPart();
        void loadMDPA(const std::string& rMDPAFilePath);
        void loadSettingsParameters(const std::string& rJSONFilePath);
        void initDofs();
        void initProperties();
        void initSolver();
        void solve();
        Kratos::ModelPart& GetMainModelPart();
        Kratos::ModelPart& GetSkinModelPart();
        Kratos::Parameters GetSettings();

    private:
        Kratos::Kernel mKernel;
        Kratos::KratosStructuralMechanicsApplication mApplication;
        std::string mModelpartName;
        Kratos::Model mModel;
        Kratos::Parameters mSettingsParameters;
        ResidualBasedNewtonRaphsonStrategyType::Pointer pmStrategy;

        Kratos::Parameters GetDefaultSettings();
    };
}
