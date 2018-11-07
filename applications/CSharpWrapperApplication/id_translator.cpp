#include "id_translator.h"

using namespace CSharpKratosWrapper;

//nodes - sorted vector with kratos IDs of surface nodes.
void IdTranslator::init(std::vector<int>& nodes) {
	int nodesSize = nodes.size();
	pmKratosIds = new int[nodesSize];

	for (int i = 0; i < nodesSize; i++) {
		pmKratosIds[i] = nodes.at(i);
	}
	pmUnityIds = new int[pmKratosIds[nodesSize - 1]];
	for (int i = 0; i < nodesSize; i++) {
		pmUnityIds[pmKratosIds[i]] = i;
	}
}

int IdTranslator::getUnityId(int kratosId) {
	return pmUnityIds[kratosId];
}

int IdTranslator::getKratosId(int unityId) {
	return pmKratosIds[unityId];
}