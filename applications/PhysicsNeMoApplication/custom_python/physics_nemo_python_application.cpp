//  ______  _                _            ______        ______                           _ _                 _
// (_____ \| |              (_)          |  ___ \      |  ___ \         /\              | (_)           _   (_)
//  _____) ) | _  _   _  ___ _  ____  ___| |   | | ____| | _ | | ___   /  \  ____  ____ | |_  ____ ____| |_  _  ___  ____
// |  ____/| || \| | | |/___) |/ ___)/___) |   | |/ _  ) || || |/ _ \ / /\ \|  _ \|  _ \| | |/ ___) _  |  _)| |/ _ \|  _ \
// | |     | | | | |_| |___ | ( (___|___ | |   | ( (/ /| || || | |_| | |__| | | | | | | | | ( (__( ( | | |__| | |_| | | | |
// |_|     |_| |_|\__  (___/|_|\____|___/|_|   |_|\____)_||_||_|\___/|______| ||_/| ||_/|_|_|\____)_||_|\___)_|\___/|_| |_|
//               (____/                                                     |_|   |_|
//
//  License:         BSD License
//                   license: PhysicsNeMoApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/define_python.h"
#include "physics_nemo_application.h"
#include "physics_nemo_application_variables.h"

namespace Kratos::Python {

PYBIND11_MODULE(KratosPhysicsNeMoApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosPhysicsNeMoApplication,
        KratosPhysicsNeMoApplication::Pointer,
        KratosApplication>(m, "KratosPhysicsNeMoApplication")
        .def(py::init<>())
        ;
}

} // namespace Kratos::Python

#endif // KRATOS_PYTHON defined
