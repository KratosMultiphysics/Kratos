


#include <vector>
#include <algorithm>
#include <iostream>
#include "mesh_converter.h"
#include "containers/array_1d.h"
#include "vector3.h"

using namespace CSharpKratosWrapper;
using namespace std;

struct element {
	int nodes[4];
	Kratos::shared_ptr<Kratos::Element> pKratosElement;
};

struct node {
	std::vector<element> elements;
};

void convert(std::vector<element>& elements, std::vector<node>& nodes, std::vector<Kratos::shared_ptr<Kratos::Element>>& kratosElements) {
	for (auto& kratosElement : kratosElements) {
		Kratos::GeometricalObject::GeometryType* pElementGeometry = &kratosElement->GetGeometry();
		struct element e = {
			{ static_cast<int>(pElementGeometry->GetPoint(0).GetId()),
				static_cast<int>(pElementGeometry->GetPoint(1).GetId()),
				static_cast<int>(pElementGeometry->GetPoint(2).GetId()),
				static_cast<int>(pElementGeometry->GetPoint(3).GetId())
			}, kratosElement
		};
		std::sort(std::begin(e.nodes), std::end(e.nodes));
		for (int i = 0; i < 4; i++) {
			if (nodes.size() <= e.nodes[i]) {
				nodes.resize(e.nodes[i] + 1);
			}
			nodes[e.nodes[i]].elements.push_back(e);
		}
		elements.push_back(e);
	}
}

bool checkContains(element& e, int(&face)[4]) {
	
	int f = 0;
	for (int i = 0; i < 4; i++) {
		if (e.nodes[i] == face[f]) {
			f++;
			if (f >= 3) {
				return true;
			}
		}
		if (i == 1 && f == 0) return false;
	}
	return false;
}

//Ensure, that node order is clockwise
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

	//Given tertrahedra ABCD, check sign([ABxBC]*AD)

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

	if (dot > 0) {
		double tmp = face.nodes[0];
		face.nodes[0] = face.nodes[1];
		face.nodes[1] = tmp;
	}
}

void process(std::vector<element>& elements, std::vector<node>& nodes, std::vector<face>& result) {
	for (auto& e : elements) {
		
		//Four faces of given tetrahedra
		int faces[4][4] = {
			{ e.nodes[0], e.nodes[1], e.nodes[2], e.nodes[3] },
		{ e.nodes[0], e.nodes[1], e.nodes[3], e.nodes[2] },
		{ e.nodes[0], e.nodes[2], e.nodes[3], e.nodes[1] },
		{ e.nodes[1], e.nodes[2], e.nodes[3], e.nodes[0] },
		};


		for (int f = 0; f < 4; f++) {

			//For each element containing first node of a face, check if that element contains given face
			bool contains = false;
			for (auto& toCheck : nodes[faces[f][0]].elements) {
				if (toCheck.pKratosElement != e.pKratosElement)
					contains = checkContains(toCheck, faces[f]);
				if (contains) break;
			}
			if (!contains) {
				face resultFace = { { faces[f][0], faces[f][1], faces[f][2], faces[f][3] } };
				fixFace(resultFace, e.pKratosElement);
				result.push_back(resultFace);
			}
		}
	}
}
int findMaxNode(std::vector<element>& elements) {
	int max = -1;
	for (auto& element : elements)
		if (element.nodes[3] > max) max = element.nodes[3];
	return max;
}

//Create sorted vector of surface node IDs
void extractNodes(std::vector<face>& faces, std::vector<int>& nodes, int maxNode) {
	bool* nodeFlags = new bool[maxNode + 1];
	for (int i = 0; i < maxNode + 1; i++)nodeFlags[i] = false;

	//mark nodes that are used at least once
	for (auto& face : faces) {
		for (int i = 0; i < 3; i++)
			nodeFlags[face.nodes[i]] = true;
	}

	//save marked nodes
	for (int i = 0; i < maxNode + 1; i++)
		if (nodeFlags[i])
			nodes.push_back(i);

	delete[] nodeFlags;
}

int findNode(int toFind, std::vector<int>& nodes) {
	int start = 0;
	int end = nodes.size() - 1;
	while (start <= end)
	{
		int mid = start + ((end - start) >> 1);

		if (nodes[mid] == toFind)
			return mid;
		else if (nodes[mid] > toFind)
			end = mid - 1;
		else if (nodes[mid] < toFind)
			start = mid + 1;
	}
	return -1;
}

//Translate node IDs into Unity format
void translateFaces(std::vector<face>& faces, std::vector<int>& nodes) {
	for (auto& f : faces) {
		for (int i = 0; i < 3; i++)f.nodes[i] = findNode(f.nodes[i], nodes);
	}
}

void MeshConverter::ProcessMesh(std::vector<Kratos::shared_ptr<Kratos::Element>>& kratosElements) {
	std::vector<element> elements;
	std::vector<node> nodes;
	convert(elements, nodes, kratosElements);
	process(elements, nodes, mFaces);
	int maxNode = findMaxNode(elements);
	extractNodes(mFaces, mNodes, maxNode);
	translateFaces(mFaces, mNodes);
}



std::vector<face>& MeshConverter::GetFaces() {
	return mFaces;
}

std::vector<int>& MeshConverter::GetNodes() {
	return mNodes;
}
