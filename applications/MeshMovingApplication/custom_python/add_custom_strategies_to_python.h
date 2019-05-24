//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED)
#define KRATOS_STRATEGIES_PYTHON_H_INCLUDED

// System includes
#include "includes/define_python.h"

namespace Kratos {
namespace Python {

void AddCustomStrategiesToPython(pybind11::module& m);

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined
