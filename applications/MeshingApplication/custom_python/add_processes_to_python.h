// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos::Python
{

void  AddProcessesToPython(pybind11::module& m);

}  // namespace Kratos::Python.
