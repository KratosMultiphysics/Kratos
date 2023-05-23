//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Fabian Meister
//

#if !defined(KRATOS_ADD_INPUT_OUTPUT_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_INPUT_OUTPUT_TO_PYTHON_H_INCLUDED

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"

// ==============================================================================

namespace Kratos
{
        namespace Python
        {

                void AddCustomInputOutputToPython(pybind11::module &m);

        } // namespace Python.
} // namespace Kratos.

#endif // KRATOS_ADD_CONTROLS_TO_PYTHON_H_INCLUDED  defined
