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
#include "custom_includes/vector3.h"

namespace CSharpKratosWrapper {
    Vector3::Vector3(double xVal, double yVal, double zVal)
        : x(xVal), y(yVal), z(zVal) {}

    Vector3::Vector3()
        : x(0), y(0), z(0) {}

    Vector3* Vector3::add(double otherX, double otherY, double otherZ) {
        x += otherX;
        y += otherY;
        z += otherZ;
        return this;

    }

    Vector3* Vector3::add(Vector3& vec) {
        return add(vec.x, vec.y, vec.z);
    }

    Vector3* Vector3::sub(double otherX, double otherY, double otherZ) {
        x -= otherX;
        y -= otherY;
        z -= otherZ;
        return this;
    }

    Vector3* Vector3::sub(Vector3& vec) {
        return sub(vec.x, vec.y, vec.z);
    }



    Vector3* Vector3::cross(Vector3& vec) {
        return new Vector3(
            y*vec.z - z * vec.y,
            z*vec.x - x * vec.z,
            x*vec.y - y * vec.x
        );
    }

    double Vector3::dot(Vector3& vec) {
        return x * vec.x + y * vec.y + z * vec.z;
    }
}
