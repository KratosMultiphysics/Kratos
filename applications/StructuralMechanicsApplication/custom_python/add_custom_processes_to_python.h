// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#pragma once

// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{
	
void  AddCustomProcessesToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.
