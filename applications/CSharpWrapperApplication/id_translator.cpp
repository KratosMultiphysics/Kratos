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
