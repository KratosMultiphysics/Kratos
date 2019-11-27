//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:    Hubert Balcerzak
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_includes/id_translator.h"

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
