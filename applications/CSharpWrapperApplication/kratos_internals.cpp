#include "kratos_internals.h"
#include "utilities/variable_utils.h"
#include "includes/model_part.h"
#include "containers/model.h"


#define SKIN_SUBMODEL_PART_NAME "main_model_part.skin_model_part"
#define MAIN_MODEL_PART_NAME "main_model_part"

using namespace CSharpKratosWrapper;

void KratosInternals::initInternals() {
	mApplication.Register();
	mKernel.Initialize();
}

void KratosInternals::loadMDPA(const std::string &mdpaPath) {
	pmModel = Kratos::make_shared<Kratos::Model>();
	Kratos::ModelPart &rMainModelPart = pmModel->CreateModelPart(MAIN_MODEL_PART_NAME, 1);
	rMainModelPart.AddNodalSolutionStepVariable(Kratos::DISPLACEMENT);
	rMainModelPart.AddNodalSolutionStepVariable(Kratos::REACTION);
	rMainModelPart.AddNodalSolutionStepVariable(Kratos::VOLUME_ACCELERATION);

	Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
	pFile->open(mdpaPath, std::fstream::in);
	Kratos::ModelPartIO(pFile).ReadModelPart(rMainModelPart);
	pFile->close();

	Kratos::ConstitutiveLaw::Pointer pCl = Kratos::make_shared<Kratos::HyperElasticIsotropicNeoHookean3D>();
	rMainModelPart.GetProperties(0).SetValue(Kratos::CONSTITUTIVE_LAW, pCl);

	rMainModelPart.SetBufferSize(2);


	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_X, Kratos::REACTION_X, rMainModelPart);
	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Y, Kratos::REACTION_Y, rMainModelPart);
	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Z, Kratos::REACTION_Z, rMainModelPart);


	rMainModelPart.CreateSubModelPart("skin_model_part");

}

void KratosInternals::initSolver() {
	SkylineLUFactorizationSolverType::Pointer pSolver = Kratos::make_shared<SkylineLUFactorizationSolverType>();
	ResidualBasedEliminationBuilderAndSolverType::Pointer pBuilderAndSolver = Kratos::make_shared<ResidualBasedEliminationBuilderAndSolverType>(
			pSolver);
	ResidualBasedIncrementalUpdateStaticSchemeType::Pointer pScheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticSchemeType>();
	ResidualCriteriaType::Pointer pConvergenceCriterion = Kratos::make_shared<ResidualCriteriaType>(1e-14, 1e-20);
	pConvergenceCriterion->SetEchoLevel(0);

	int maxIters = 20;
	bool computerReactions = true;
	bool reformStepDofs = true;
	bool moveMeshFlag = true;

	pmStrategy = Kratos::make_shared<ResidualBasedNewtonRaphsonStrategyType>(
			pmModel->GetModelPart(MAIN_MODEL_PART_NAME),
			pScheme,
			pSolver,
			pConvergenceCriterion,
			pBuilderAndSolver,
			maxIters,
			computerReactions,
			reformStepDofs,
			moveMeshFlag);

	pmStrategy->SetEchoLevel(0);
	pmStrategy->Check();
}

void KratosInternals::solve() {
	pmStrategy->Solve();
}

Kratos::ModelPart &KratosInternals::rGetMainModelPart() {
	return pmModel->GetModelPart(MAIN_MODEL_PART_NAME);
}

Kratos::ModelPart &KratosInternals::rGetSkinModelPart() {
	return pmModel->GetModelPart(SKIN_SUBMODEL_PART_NAME);
}