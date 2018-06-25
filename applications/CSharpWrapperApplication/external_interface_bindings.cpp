#define EXPORT __declspec(dllexport)
#include <stdlib.h>
#include "external_interface.h"
#include <iostream>

using namespace std;
using namespace CSharpKratosWrapper;

extern "C" {

	EXPORT void __stdcall Init(char* path) {
		Interface::init(path);
	}

	EXPORT float* __stdcall GetXCoordinates() {
		return Interface::getXCoordinates();
	}

	EXPORT float* __stdcall GetYCoordinates() {
		return Interface::getYCoordinates();
	}

	EXPORT float* __stdcall GetZCoordinates() {
		return Interface::getZCoordinates();
	}

	EXPORT int __stdcall GetNodesCount() {
		return Interface::getNodesCount();
	}

	EXPORT int* __stdcall GetTriangles() {
		return Interface::getTriangles();
	}

	EXPORT int __stdcall GetTrianglesCount() {
		return Interface::getTrianglesCount();
	}

	EXPORT void __stdcall UpdateNodePos(int nodeId, float x, float y, float z) {
		Interface::updateNodePos(nodeId, x, y, z);
	}

	EXPORT void __stdcall Calculate() {
		Interface::calculate();
	}
}