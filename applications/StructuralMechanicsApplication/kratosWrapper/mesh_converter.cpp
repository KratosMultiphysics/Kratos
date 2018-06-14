#include <vector>
#include <algorithm>
#include <iostream>
#include "mesh_converter.h"
#include "containers/array_1d.h"
#include "vector3.h"

using namespace KratosWrapper;

struct element {
	int nodes[4];
	Kratos::shared_ptr<Kratos::Element> pKratosElement;
};

void convert(std::vector<element>& converted, std::vector<Kratos::shared_ptr<Kratos::Element>>& input) {
	for (auto &kratosElement : input) {
		Kratos::GeometricalObject::GeometryType* pElementGeometry = &kratosElement->GetGeometry();
		struct element e = { { static_cast<int>(pElementGeometry->GetPoint(0).GetId()),
			static_cast<int>(pElementGeometry->GetPoint(1).GetId()),
			static_cast<int>(pElementGeometry->GetPoint(2).GetId()),
			static_cast<int>(pElementGeometry->GetPoint(3).GetId())
			}, kratosElement };
		converted.push_back(e);
	}
}

bool elementSorter(element const& e1, element const& e2) {
	for (int i = 0; i < 4; i++) {
		if (e1.nodes[i] != e2.nodes[i]) return e1.nodes[i] < e2.nodes[i];
	}
	return false;
}

void sortElements(std::vector<element>& elements) {
	for (auto &e : elements) {
		std::sort(std::begin(e.nodes), std::end(e.nodes));
	}
	std::sort(elements.begin(), elements.end(), elementSorter);
}

bool checkContains(int(&el)[4], int(&face)[4], bool* faceWasUsed) {
	int f = 0;
	int nface = -1;
	for (int i = 0; i < 4; i++) {
		if (el[i] == face[f]) {
			if (nface == -1) {
				if (f == 0 && i == 1)nface = 3;
				else if (f == 1 && i == 2)nface = 2;
				else if (f == 2 && i == 3)nface = 1;
				else if (f == 2 && i == 2)nface = 0;
			}
			f++;
			if (f >= 3) {
				faceWasUsed[nface] = true;
				return true;
			}
		}
		if (i == 1 && f == 0) return false;
	}
	return false;
}



void fixFace(face& face, Kratos::shared_ptr < Kratos::Element > kratosElement) {
	Kratos::array_1d<double, 3> points[4];

	//assign coordinates
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (face.nodes[i] == kratosElement->pGetGeometry()->pGetPoint(j)->Id()) {
				points[i] = kratosElement->pGetGeometry()->pGetPoint(j)->Coordinates();
				break;
			}
		}

	}
	Vector3* AB = new Vector3(points[1][0], points[1][1], points[1][2]);
	AB->sub(points[0][0], points[0][1], points[0][2]);

	Vector3* BC = new Vector3(points[2][0], points[2][1], points[2][2]);
	BC->sub(points[1][0], points[1][1], points[1][2]);
	
	Vector3* AD = new Vector3(points[3][0], points[3][1], points[3][2]);
	AD->sub(points[0][0], points[0][1], points[0][2]);

	Vector3* normal = AB->cross(*BC);

	double dot = normal->dot(*AD);
	delete AB;
	delete BC;
	delete AD;
	delete normal;

	if (dot < 0) {
		double tmp = face.nodes[0];
		face.nodes[0] = face.nodes[1];
		face.nodes[1] = tmp;
	}
}

void process(std::vector<element>& elements, std::vector<face>& result) {
	bool** wasUsed = new bool*[elements.size()];

	for (int i = 0; i < elements.size(); i++) {
		wasUsed[i] = new bool[4];
		for (int j = 0; j < 4; j++)wasUsed[i][j] = false;
	}

	for (int i = 0; i < elements.size(); i++) {
		element& e = elements[i];
		int faces[4][4] = {
			{ e.nodes[0], e.nodes[1], e.nodes[2], e.nodes[3] },
		{ e.nodes[0], e.nodes[1], e.nodes[3], e.nodes[2] },
		{ e.nodes[0], e.nodes[2], e.nodes[3], e.nodes[1] },
		{ e.nodes[1], e.nodes[2], e.nodes[3], e.nodes[0] },
		};

		for (int f = 0; f < 4; f++) {
			if (!wasUsed[i][f]) {
				int current = i + 1;
				bool faceFound = false;
				while (!faceFound && current < elements.size() && elements[current].nodes[0] <= faces[f][0]) {
					faceFound = checkContains(elements[current].nodes, faces[f], wasUsed[current]);
					current++;
				}

				if (!faceFound) {
					face resultFace = { { faces[f][0], faces[f][1], faces[f][2], faces[f][3] } };
					fixFace(resultFace, e.pKratosElement);
					result.push_back(resultFace);
				}
			}
		}
	}

	for (int i = 0; i < elements.size(); i++) {
		delete wasUsed[i];
	}
	delete wasUsed;
}

int findMaxNode(std::vector<element> elements) {
	int max = -1;
	for (auto& element : elements)
		if (element.nodes[3] > max) max = element.nodes[3];
	return max;
}

void extractNodes(std::vector<face>& faces, std::vector<int>& nodes, int maxNode) {
	bool* nodeFlags = new bool[maxNode+1];
	for (int i = 0; i < maxNode+1; i++)nodeFlags[i] = false;
	
	//mark nodes that are used at least once
	for (auto& face : faces) {
		for (int i = 0; i < 3; i++)
			nodeFlags[face.nodes[i]] = true;
	}

	//save marked nodes
	for (int i = 0; i < maxNode + 1; i++)
		if (nodeFlags[i])
			nodes.push_back(i);
	
	delete nodeFlags;
}

int findNode(int toFind, std::vector<int>& nodes) {
	int mid;
	int start = 0;
	int end = nodes.size()-1;
	while (start <= end)
	{
		mid = start + ((end - start) >> 1);

		if (nodes[mid] == toFind)
			return mid;
		else if (nodes[mid] > toFind)
			end = mid - 1;
		else if (nodes[mid] < toFind)
			start = mid + 1;
	}
}

void translateFaces(std::vector<face>& faces, std::vector<int>& nodes) {
	for (auto& f : faces) {
		for (int i = 0; i < 3; i++)f.nodes[i] = findNode(f.nodes[i], nodes);
	}
}

void printFaces(std::vector<face>& faces) {
	for (auto& f : faces) {
		std::cout << f.nodes[0] << " " << f.nodes[1] << " " << f.nodes[2] << std::endl;
	}
	std::cout << "---------------------------" << std::endl;

}

void MeshConverter::ProcessMesh(std::vector<Kratos::shared_ptr<Kratos::Element>>& kratosElements) {
	std::vector<element> elements;
	convert(elements, kratosElements);
	sortElements(elements);
	process(elements, mFaces);
	int maxNode = findMaxNode(elements);
	extractNodes(mFaces, mNodes, maxNode);
	//printFaces(mFaces);
	translateFaces(mFaces, mNodes);
	//printFaces(mFaces);
}



std::vector<face>& MeshConverter::GetFaces() {
	return mFaces;
}

std::vector<int>& MeshConverter::GetNodes() {
	return mNodes;
}

