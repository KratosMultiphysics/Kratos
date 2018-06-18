#include "includes/model_part_io.h"
#include "includes/variables.h"
#include "external_interface.h"

using namespace KratosWrapper;


int* Interface::pmKratosNodeIds;
int* Interface::pmTriangles;
int Interface::mTrianglesCount;
int Interface::mNodesCount;
float* Interface::pmXCoordinates;
float* Interface::pmYCoordinates;
float* Interface::pmZCoordinates;
Kratos::Kernel Interface::mKernel;
Kratos::KratosStructuralMechanicsApplication Interface::mApplication;
Kratos::ModelPart Interface::mMainModelPart;

void Interface::initInternals() {
	mApplication.Register();
	mKernel.Initialize();
	mMainModelPart.AddNodalSolutionStepVariable(Kratos::DISPLACEMENT);
}

void Interface::loadMDPAFile(char* mdpaPath){
	Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
	pFile->open(mdpaPath, std::fstream::in);

	Kratos::ModelPartIO(pFile).ReadModelPart(mMainModelPart);
	pFile->close();
}

void Interface::saveTriangles(MeshConverter& meshConverter) {
	std::vector<face> faces = meshConverter.GetFaces();

	mTrianglesCount = faces.size();
	pmTriangles = new int[mTrianglesCount * 3];
	for (int i = 0; i < mTrianglesCount; i++) {
		for (int j = 0; j < 3; j++) {
			pmTriangles[3 * i + j] = faces.at(i).nodes[j];
		}
	}
}

void Interface::saveNodes(MeshConverter& meshConverter) {
	std::vector<int> nodes = meshConverter.GetNodes();

	//save nodes
	mNodesCount = nodes.size();
	pmKratosNodeIds = new int[mNodesCount];
	for (int i = 0; i < mNodesCount; i++) {
		pmKratosNodeIds[i] = nodes.at(i);
	}

	pmXCoordinates = new float[mNodesCount];
	pmYCoordinates = new float[mNodesCount];
	pmZCoordinates = new float[mNodesCount];

	//retrieve node coordinates
	for (int i = 0; i < mNodesCount; i++) {
		Kratos::shared_ptr<Kratos::Node<3Ui64>> currentNode = mMainModelPart.pGetNode(pmKratosNodeIds[i]);
		pmXCoordinates[i] = currentNode->X();
		pmYCoordinates[i] = currentNode->Y();
		pmZCoordinates[i] = currentNode->Z();

	}
}

void Interface::init(char* mdpaPath) {
	initInternals();
	loadMDPAFile(mdpaPath);

	MeshConverter meshConverter;
	meshConverter.ProcessMesh(mMainModelPart.ElementsArray());

	saveTriangles(meshConverter);
	saveNodes(meshConverter);
}

void Interface::updateNodePos(int nodeId, float x, float y, float z) {
	Kratos::NodeType& node =  mMainModelPart.GetNode(pmKratosNodeIds[nodeId]);
	Kratos::array_1d<double, 3>& displacement = node.FastGetSolutionStepValue(Kratos::DISPLACEMENT);
	displacement[0] = x - node.X();
	displacement[1] = y - node.Y();
	displacement[2] = z - node.Z();
}

void Interface::calculate() {
	//haha no
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
