#include "id_translator.h"

using namespace CSharpKratosWrapper;

//nodes - sorted vector with kratos IDs of surface nodes.
void IdTranslator::init(std::vector<int> &nodes) {
    int nodesSize = nodes.size();
    pmKratosIds = new int[nodesSize];

    for (int i = 0; i < nodesSize; i++) {
        pmKratosIds[i] = nodes.at(i);
    }
    pmSurfaceIds = new int[pmKratosIds[nodesSize - 1]];
    for (int i = 0; i < pmKratosIds[nodesSize - 1]; i++) pmSurfaceIds[i] = -1;

    for (int i = 0; i < nodesSize; i++) {
        pmSurfaceIds[pmKratosIds[i]] = i;
    }
}

int IdTranslator::getSurfaceId(int kratosId) {
    return pmSurfaceIds[kratosId];
}

int IdTranslator::getKratosId(int unityId) {
    return pmKratosIds[unityId];
}

bool IdTranslator::hasKratosId(int surfaceId) {
    return surfaceId < mKratosIdsSize;
}

bool IdTranslator::hasSurfaceId(int kratosId) {
    return mSurfaceIdsSize < kratosId && pmSurfaceIds[kratosId] != -1;
}

int IdTranslator::safeGetKratosId(int surfaceId) {
    if (surfaceId >= mKratosIdsSize) return -1;
    return pmKratosIds[surfaceId];
}

int IdTranslator::safeGetSurfaceId(int kratosId) {
    if (kratosId >= mSurfaceIdsSize) return -1;
    return pmSurfaceIds[kratosId];
}
