//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#ifdef KRATOS_PYTHON
#include <pybind11/pybind11.h>

#include "custom_python/add_trilinos_schemes_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosPoromechanicsTrilinosExtension, m)
{
    AddTrilinosSchemesToPython(m);
}

}  // namespace Python.
}  // namespace Kratos.

#endif
