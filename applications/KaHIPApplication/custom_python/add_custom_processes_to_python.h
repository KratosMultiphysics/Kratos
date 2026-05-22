//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python {

/**
 * @brief Register KaHIPApplication process classes with the Python module.
 * @details Exposes:
 *          - KaHIPDivideHeterogeneousInputProcess
 *          - KaHIPDivideHeterogeneousInputInMemoryProcess
 * @param m  The pybind11 module object
 */
void AddCustomProcessesToPython(pybind11::module& m);

} // namespace Kratos::Python
