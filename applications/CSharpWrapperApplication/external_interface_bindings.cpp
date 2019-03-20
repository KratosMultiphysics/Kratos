#define EXPORT __declspec(dllexport)
#include <stdlib.h>
#include "external_interface.h"
#include <iostream>

using namespace std;
using namespace CSharpKratosWrapper;

extern "C" {

	EXPORT void __stdcall Init(char* path) {
		CSharpInterface::init(path);
	}

	EXPORT float* __stdcall GetXCoordinates() {
		return CSharpInterface::getXCoordinates();
	}

	EXPORT float* __stdcall GetYCoordinates() {
		return CSharpInterface::getYCoordinates();
	}

	EXPORT float* __stdcall GetZCoordinates() {
		return CSharpInterface::getZCoordinates();
	}

	EXPORT int __stdcall GetNodesCount() {
		return CSharpInterface::getNodesCount();
	}

	EXPORT int* __stdcall GetTriangles() {
		return CSharpInterface::getTriangles();
	}

	EXPORT int __stdcall GetTrianglesCount() {
		return CSharpInterface::getTrianglesCount();
	}

	EXPORT void __stdcall UpdateNodePos(int nodeId, float x, float y, float z) {
		CSharpInterface::updateNodePos(nodeId, x, y, z);
	}

	EXPORT void __stdcall Calculate() {
		CSharpInterface::calculate();
	}
}