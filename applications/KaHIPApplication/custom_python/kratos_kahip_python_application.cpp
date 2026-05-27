//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "kahip_application.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_modeler_to_python.h"

namespace Kratos::Python {

PYBIND11_MODULE(KratosKaHIPApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosKaHIPApplication,
               KratosKaHIPApplication::Pointer,
               KratosApplication>(m, "KratosKaHIPApplication")
        .def(py::init<>())
        ;

    AddCustomProcessesToPython(m);
    AddCustomModelerToPython(m);
}

} // namespace Kratos::Python

#endif // KRATOS_PYTHON defined
