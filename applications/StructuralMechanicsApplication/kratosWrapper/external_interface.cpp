#include "includes/variables.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "utilities/variable_utils.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

#include "external_interface.h"

#include "includes/constitutive_law.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

using namespace KratosWrapper;


int* Interface::pmKratosNodeIds;
int* Interface::pmUnityNodeIds;
int* Interface::pmTriangles;
int Interface::mTrianglesCount;
int Interface::mNodesCount;
float* Interface::pmXCoordinates;
float* Interface::pmYCoordinates;
float* Interface::pmZCoordinates;
Kratos::Kernel Interface::mKernel;
Kratos::KratosStructuralMechanicsApplication Interface::mApplication;
Kratos::ModelPart Interface::mMainModelPart;
std::vector<Kratos::NodeType::Pointer> Interface::mFixedNodes;

void Interface::initInternals() {
	mApplication.Register();
	mKernel.Initialize();
}

void Interface::loadMDPAFile(char* mdpaPath){
	mMainModelPart.AddNodalSolutionStepVariable(Kratos::DISPLACEMENT);
	mMainModelPart.AddNodalSolutionStepVariable(Kratos::REACTION);
	mMainModelPart.AddNodalSolutionStepVariable(Kratos::VOLUME_ACCELERATION);
	
	Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
	pFile->open(mdpaPath, std::fstream::in);
	Kratos::ModelPartIO(pFile).ReadModelPart(mMainModelPart);
	pFile->close();

	//LOAD THE CONSTITUTIVE LAW
	Kratos::ConstitutiveLaw::Pointer pCl = Kratos::make_shared<Kratos::HyperElasticIsotropicNeoHookean3D>();
	mMainModelPart.GetProperties(0).SetValue(Kratos::CONSTITUTIVE_LAW, pCl);
	mMainModelPart.SetBufferSize(2);

}

void Interface::saveTriangles(MeshConverter& meshConverter) {

	Kratos::ModelPart::Pointer skinModelPart = mMainModelPart.CreateSubModelPart("skin_model_part");
	int lastId = mMainModelPart.Elements().back().Id();
	
	std::vector<face> faces = meshConverter.GetFaces();

	mTrianglesCount = faces.size();
	pmTriangles = new int[mTrianglesCount * 3];
	
	for (int i = 0; i < mTrianglesCount; i++) {
		std::vector<Kratos::IndexType> nodes;
		for (int j = 0; j < 3; j++) {
			pmTriangles[3 * i + j] = faces.at(i).nodes[j];
			nodes.push_back(pmKratosNodeIds[faces.at(i).nodes[j]]);
			skinModelPart->AddNode(mMainModelPart.pGetNode(pmKratosNodeIds[faces.at(i).nodes[j]]));
		}
		lastId++;
		skinModelPart->CreateNewCondition("SurfaceCondition3D3N", lastId, nodes, mMainModelPart.pGetProperties(0));
	}
}

void Interface::saveNodes(MeshConverter& meshConverter) {
	std::vector<int> nodes = meshConverter.GetNodes();
	mNodesCount = nodes.size();
	pmKratosNodeIds = new int[mNodesCount];

	for (int i = 0; i < mNodesCount; i++) {
		pmKratosNodeIds[i] = nodes.at(i);
	}
	pmUnityNodeIds = new int[pmKratosNodeIds[mNodesCount - 1]];
	for (int i = 0; i < mNodesCount; i++) {
		pmUnityNodeIds[pmKratosNodeIds[i]] = i;
	}

	pmXCoordinates = new float[mNodesCount];
	pmYCoordinates = new float[mNodesCount];
	pmZCoordinates = new float[mNodesCount];
}


void Interface::retrieveNodesPos() {
	Kratos::ModelPart::Pointer skin_part = mMainModelPart.pGetSubModelPart("skin_model_part");

#pragma omp parallel for
	for (int i = 0; i<skin_part->Nodes().size(); ++i)
	{
		auto currentNode = skin_part->NodesBegin() + i;
		int currentNodeUnityId = pmUnityNodeIds[currentNode->Id()];
		pmXCoordinates[currentNodeUnityId] = currentNode->X();
		pmYCoordinates[currentNodeUnityId] = currentNode->Y();
		pmZCoordinates[currentNodeUnityId] = currentNode->Z();
	}
}

void Interface::init(char* mdpaPath) {

	initInternals();
	loadMDPAFile(mdpaPath);

	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_X, Kratos::REACTION_X, mMainModelPart);
	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Y, Kratos::REACTION_Y, mMainModelPart);
	Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Z, Kratos::REACTION_Z, mMainModelPart);

	MeshConverter meshConverter;
	meshConverter.ProcessMesh(mMainModelPart.ElementsArray());

	saveNodes(meshConverter);
	saveTriangles(meshConverter);
	retrieveNodesPos();
}

void Interface::updateNodePos(int nodeId, float x, float y, float z) {
	Kratos::NodeType::Pointer node =  mMainModelPart.pGetNode(pmKratosNodeIds[nodeId]);
	node->Fix(Kratos::DISPLACEMENT_X);
	node->Fix(Kratos::DISPLACEMENT_Y);
	node->Fix(Kratos::DISPLACEMENT_Z);
	Kratos::array_1d<double, 3>& displacement = node->FastGetSolutionStepValue(Kratos::DISPLACEMENT);

	displacement[0] = x - node->X0();
	displacement[1] = y - node->Y0();
	displacement[2] = z - node->Z0();
	mFixedNodes.push_back(node);
}

void Interface::freeNodes() {
	for (auto& node : mFixedNodes) {
		node->Free(Kratos::DISPLACEMENT_X);
		node->Free(Kratos::DISPLACEMENT_Y);
		node->Free(Kratos::DISPLACEMENT_Z);
	} 
	mFixedNodes.clear();
}


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


void Interface::calculate() {
	SkylineLUFactorizationSolverType::Pointer pSolver = Kratos::make_shared< SkylineLUFactorizationSolverType>();
	ResidualBasedEliminationBuilderAndSolverType::Pointer pBuilderAndSolver = Kratos::make_shared<ResidualBasedEliminationBuilderAndSolverType>(pSolver);
	ResidualBasedIncrementalUpdateStaticSchemeType::Pointer pScheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticSchemeType>();
	ResidualCriteriaType::Pointer pConvergenceCriterion = Kratos::make_shared<ResidualCriteriaType>(1e-14, 1e-20);
	pConvergenceCriterion->SetEchoLevel(0);


	//mMainModelPart.Nodes()[pmKratosNodeIds[1]].Fix(Kratos::DISPLACEMENT_X);
	//mMainModelPart.Nodes()[pmKratosNodeIds[1]].FastGetSolutionStepValue(Kratos::DISPLACEMENT_X, 0) = 0.1;
	
	int maxIters = 20;
	bool computerReactions = true;
	bool reformStepDofs = true;
	bool moveMeshFlag = true;

	ResidualBasedNewtonRaphsonStrategyType::Pointer pStrategy = Kratos::make_shared < ResidualBasedNewtonRaphsonStrategyType >(
		mMainModelPart,
		pScheme,
		pSolver,
		pConvergenceCriterion,
		pBuilderAndSolver,
		maxIters,
		computerReactions,
		reformStepDofs,
		moveMeshFlag);

	pStrategy->SetEchoLevel(0);

	pStrategy->Check();
	pStrategy->Solve();

	retrieveNodesPos();
	freeNodes();
}

float* Interface::getXCoordinates() {
	return pmXCoordinates;
}

float* Interface::getYCoordinates() {
	return pmYCoordinates;
}

float* Interface::getZCoordinates() {
	return pmZCoordinates;
}

int Interface::getNodesCount() {
	return mNodesCount;
}

int* Interface::getTriangles() {
	return pmTriangles;
}

int Interface::getTrianglesCount() {
	return mTrianglesCount;
}
