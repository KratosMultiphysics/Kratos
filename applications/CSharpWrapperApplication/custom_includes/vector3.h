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

#if !defined(CSHARP_WRAPPER_APPLICATION_VECTOR3_H_INCLUDED )
#define  CSHARP_WRAPPER_APPLICATION_VECTOR3_H_INCLUDED

// System includes

// External includes

// Project includes

namespace CSharpKratosWrapper {

    class Vector3 {
    public:
        double x;
        double y;
        double z;
        Vector3(double x, double y, double z);
        Vector3();
        Vector3* cross(Vector3&);
        double dot(Vector3&);
        Vector3* add(Vector3&);
        Vector3* add(double otherX, double otherY, double otherZ);
        Vector3* sub(Vector3&);
        Vector3* sub(double otherX, double otherY, double otherZ);
    };

}

#endif	/* CSHARP_WRAPPER_APPLICATION_VECTOR3_H_INCLUDED */
