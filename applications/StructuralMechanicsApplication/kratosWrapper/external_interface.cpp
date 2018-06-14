#include "includes/model_part_io.h"

#include "mesh_converter.h"
#include "external_interface.h"

using namespace KratosWrapper;


int* Interface::pmRequiredNodes;
int* Interface::pmTriangles;
int Interface::mTrianglesCount;
int Interface::mNodesCount;
float* Interface::pmXCoordinates;
float* Interface::pmYCoordinates;
float* Interface::pmZCoordinates;
Kratos::Kernel Interface::mKernel;
Kratos::KratosStructuralMechanicsApplication Interface::mApplication;
Kratos::ModelPart Interface::mMainModelPart;

void Interface::init(char* path) {
	mApplication.Register();
	mKernel.Initialize();

	Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
	pFile->open(path, std::fstream::in);

	Kratos::ModelPartIO * modelPartIO = new Kratos::ModelPartIO(pFile);
	modelPartIO->ReadModelPart(mMainModelPart);
	delete modelPartIO;

	MeshConverter meshConverter;
	meshConverter.ProcessMesh(mMainModelPart.ElementsArray());

	std::vector<face> faces = meshConverter.GetFaces();
	std::vector<int> nodes = meshConverter.GetNodes();
	
	//save converted traingle mesh
	mTrianglesCount = faces.size();
	pmTriangles = new int[mTrianglesCount * 3];
	for (int i = 0; i < mTrianglesCount; i++) {
		for (int j = 0; j < 3; j++) {
			pmTriangles[3 * i + j] = faces.at(i).nodes[j];
		}
	}

	//save nodes
	mNodesCount = nodes.size();
	pmRequiredNodes = new int[mNodesCount];
	for (int i = 0; i < mNodesCount; i++) {
		pmRequiredNodes[i] = nodes.at(i);
	}

	pmXCoordinates = new float[mNodesCount];
	pmYCoordinates = new float[mNodesCount];
	pmZCoordinates = new float[mNodesCount];

	//retrieve node coordinates
	for (int i = 0; i < mNodesCount; i++) {
		Kratos::shared_ptr<Kratos::Node<3Ui64>> currentNode = mMainModelPart.pGetNode(pmRequiredNodes[i]);
		pmXCoordinates[i] = currentNode->X();
		pmYCoordinates[i] = currentNode->Y();
		pmZCoordinates[i] = currentNode->Z();
	}

}

void Interface::updateNodesPos(float* xCoordinates, float* yCoordinates, float* zCoordinates) {

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
