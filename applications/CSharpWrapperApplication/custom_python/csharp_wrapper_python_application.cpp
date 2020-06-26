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

#if defined(KRATOS_PYTHON)
// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

#include "csharp_wrapper_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosCSharpWrapperApplication, m)
{
    class_<KratosCSharpWrapperApplication,
    KratosCSharpWrapperApplication::Pointer,
    KratosApplication>(m, "KratosCSharpWrapperApplication")
    .def(init<>());

    // Registering variables in python
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, )
}


}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
