#include "includes/kernel.h"
#include "structural_mechanics_application.h"
#include "includes/model_part_io.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "includes/constitutive_law.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

namespace CSharpKratosWrapper {


	typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrix;
	typedef boost::numeric::ublas::vector<double> Vector;
	typedef Kratos::UblasSpace<double, CompressedMatrix, Vector> SpaceType;
	typedef boost::numeric::ublas::matrix<double> Matrix;
	typedef Kratos::UblasSpace<double, Matrix, Vector> LocalSpaceType;
	typedef Kratos::Reorderer<SpaceType, LocalSpaceType > ReordererType;
	typedef Kratos::SkylineLUFactorizationSolver<SpaceType, LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

	typedef Kratos::UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef Kratos::LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
	typedef Kratos::ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;

	typedef Kratos::ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;

	typedef Kratos::ResidualCriteria<SparseSpaceType, LocalSpaceType > ResidualCriteriaType;

	typedef Kratos::ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;




	class KratosInternals {

	public:
		void initInternals();
		void loadMDPA(std::string mdpaPath);
		void initSolver();
		void solve();
		Kratos::ModelPart::Pointer pGetMainModelPart();
		Kratos::ModelPart::Pointer pGetSkinModelPart();

	private:
		Kratos::Kernel mKernel;
		Kratos::KratosStructuralMechanicsApplication mApplication;
		Kratos::ModelPart::Pointer pmMainModelPart;
		ResidualBasedNewtonRaphsonStrategyType::Pointer pmStrategy;
	};
}