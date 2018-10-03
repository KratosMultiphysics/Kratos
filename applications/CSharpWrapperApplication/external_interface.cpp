#include "external_interface.h"

#include "includes/constitutive_law.h"

using namespace CSharpKratosWrapper;

int* CSharpInterface::pmTriangles;
int CSharpInterface::mTrianglesCount;
int CSharpInterface::mNodesCount;
float* CSharpInterface::pmXCoordinates;
float* CSharpInterface::pmYCoordinates;
float* CSharpInterface::pmZCoordinates;
KratosInternals CSharpInterface::mKratosInternals;
std::vector<Kratos::NodeType::Pointer> CSharpInterface::mFixedNodes;
IdTranslator CSharpInterface::idTranslator;

//Fetch faces from mesh converter and create skin model part
void CSharpInterface::saveTriangles(MeshConverter& meshConverter) {

	Kratos::ModelPart::Pointer pSkinModelPart = mKratosInternals.pGetSkinModelPart();
	Kratos::ModelPart::Pointer pMainModelPart = mKratosInternals.pGetMainModelPart();
	int lastId = pMainModelPart->Elements().back().Id();
	
	std::vector<face> faces = meshConverter.GetFaces();

	mTrianglesCount = faces.size();
	pmTriangles = new int[mTrianglesCount * 3];
	
	for (int i = 0; i < mTrianglesCount; i++) {
		std::vector<Kratos::IndexType> nodes;
		for (int j = 0; j < 3; j++) {
			pmTriangles[3 * i + j] = faces.at(i).nodes[j];
			nodes.push_back(idTranslator.getKratosId(faces.at(i).nodes[j]));
			pSkinModelPart->AddNode(pMainModelPart->pGetNode(idTranslator.getKratosId(faces.at(i).nodes[j])));
		}
		lastId++;
		pSkinModelPart->CreateNewCondition("SurfaceCondition3D3N", lastId, nodes, pMainModelPart->pGetProperties(0));
	}
}

//Fetch surface nodes from mesh converter and initialize ID translator
void CSharpInterface::saveNodes(MeshConverter& meshConverter) {
	std::vector<int> nodes = meshConverter.GetNodes();
	mNodesCount = nodes.size();
	
	idTranslator.init(nodes);

	pmXCoordinates = new float[mNodesCount];
	pmYCoordinates = new float[mNodesCount];
	pmZCoordinates = new float[mNodesCount];
	retrieveNodesPos();
}

//Save recalculated surface nodes positions
void CSharpInterface::retrieveNodesPos() {
	Kratos::ModelPart::Pointer skin_part = mKratosInternals.pGetSkinModelPart();

#pragma omp parallel for
	for (int i = 0; i<skin_part->Nodes().size(); ++i)
	{
		auto currentNode = skin_part->NodesBegin() + i;
		int currentNodeUnityId = idTranslator.getUnityId(currentNode->Id());
		pmXCoordinates[currentNodeUnityId] = currentNode->X();
		pmYCoordinates[currentNodeUnityId] = currentNode->Y();
		pmZCoordinates[currentNodeUnityId] = currentNode->Z();
	}
}

void CSharpInterface::freeNodes() {
	for (auto& node : mFixedNodes) {
		node->Free(Kratos::DISPLACEMENT_X);
		node->Free(Kratos::DISPLACEMENT_Y);
		node->Free(Kratos::DISPLACEMENT_Z);
	}
	mFixedNodes.clear();
}

void CSharpInterface::init(char* mdpaPath) {
	mKratosInternals.initInternals();
	mKratosInternals.loadMDPA(std::string(mdpaPath));
	mKratosInternals.initSolver();

	MeshConverter meshConverter;
	meshConverter.ProcessMesh(mKratosInternals.pGetMainModelPart()->ElementsArray());

	saveNodes(meshConverter);
	saveTriangles(meshConverter);
}

//Update DISPLACEMENT variable of a node, so that final position is as given. X0 + DISPLACEMENT_X = x
void CSharpInterface::updateNodePos(int nodeId, float x, float y, float z) {
	Kratos::NodeType::Pointer node =  mKratosInternals.pGetMainModelPart()->pGetNode(idTranslator.getKratosId(nodeId));
	node->Fix(Kratos::DISPLACEMENT_X);
	node->Fix(Kratos::DISPLACEMENT_Y);
	node->Fix(Kratos::DISPLACEMENT_Z);
	Kratos::array_1d<double, 3>& displacement = node->FastGetSolutionStepValue(Kratos::DISPLACEMENT);

	displacement[0] = x - node->X0();
	displacement[1] = y - node->Y0();
	displacement[2] = z - node->Z0();
	mFixedNodes.push_back(node);
}

void CSharpInterface::calculate() {
	mKratosInternals.solve();
	retrieveNodesPos();

	freeNodes();
	
}

float* CSharpInterface::getXCoordinates() {
	return pmXCoordinates;
}

float* CSharpInterface::getYCoordinates() {
	return pmYCoordinates;
}

float* CSharpInterface::getZCoordinates() {
	return pmZCoordinates;
}

int CSharpInterface::getNodesCount() {
	return mNodesCount;
}

int* CSharpInterface::getTriangles() {
	return pmTriangles;
}

int CSharpInterface::getTrianglesCount() {
	return mTrianglesCount;
}
