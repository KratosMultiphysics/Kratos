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
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <stdlib.h>
#include <iostream>
#define EXPORT __declspec(dllexport)

// External includes

// Project includes
#include "external_interface.h"

using namespace std;
using namespace CSharpKratosWrapper;

extern "C" {

#if defined(KRATOS_COMPILED_IN_WINDOWS)
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
#endif
}
