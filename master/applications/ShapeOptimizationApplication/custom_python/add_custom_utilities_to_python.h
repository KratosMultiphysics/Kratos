// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

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

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED  defined
