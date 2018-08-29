#include "kratos_internals.h"
#include "utilities/variable_utils.h"

#define SKIN_SUBMODEL_PART_NAME "skin_model_part"

using namespace CSharpKratosWrapper;

void KratosInternals::initInternals() {
	mApplication.Register();
	mKernel.Initialize();
}

void KratosInternals::loadMDPA(std::string mdpaPath) {
	pmMainModelPart = Kratos::make_shared<Kratos::ModelPart>();
	pmMainModelPart->AddNodalSolutionStepVariable(Kratos::DISPLACEMENT);
	pmMainModelPart->AddNodalSolutionStepVariable(Kratos::REACTION);
	pmMainModelPart->AddNodalSolutionStepVariable(Kratos::VOLUME_ACCELERATION);

	Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
	pFile->open(mdpaPath, std::fstream::in);
	Kratos::ModelPartIO(pFile).ReadModelPart(*pmMainModelPart);
	pFile->close();

	Kratos::ConstitutiveLaw::Pointer pCl = Kratos::make_shared<Kratos::HyperElasticIsotropicNeoHookean3D>();
	pmMainModelPart->GetProperties(0).SetValue(Kratos::CONSTITUTIVE_LAW, pCl);
	pmMainModelPart->SetBufferSize(2);

	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_X, Kratos::REACTION_X, *pmMainModelPart);
	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Y, Kratos::REACTION_Y, *pmMainModelPart);
	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Z, Kratos::REACTION_Z, *pmMainModelPart);

	pmMainModelPart->CreateSubModelPart(SKIN_SUBMODEL_PART_NAME);
}

void KratosInternals::initSolver() {
	SkylineLUFactorizationSolverType::Pointer pSolver = Kratos::make_shared< SkylineLUFactorizationSolverType>();
	ResidualBasedEliminationBuilderAndSolverType::Pointer pBuilderAndSolver = Kratos::make_shared<ResidualBasedEliminationBuilderAndSolverType>(pSolver);
	ResidualBasedIncrementalUpdateStaticSchemeType::Pointer pScheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticSchemeType>();
	ResidualCriteriaType::Pointer pConvergenceCriterion = Kratos::make_shared<ResidualCriteriaType>(1e-14, 1e-20);
	pConvergenceCriterion->SetEchoLevel(0);

	int maxIters = 20;
	bool computerReactions = true;
	bool reformStepDofs = true;
	bool moveMeshFlag = true;

	pmStrategy = Kratos::make_shared < ResidualBasedNewtonRaphsonStrategyType >(
		*pmMainModelPart,
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

Kratos::ModelPart::Pointer KratosInternals::pGetMainModelPart() {
	return pmMainModelPart;
}

Kratos::ModelPart::Pointer KratosInternals::pGetSkinModelPart() {
	return pmMainModelPart->pGetSubModelPart(SKIN_SUBMODEL_PART_NAME);
}