//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Ruben Zorrilla
//
//

#if !defined(KRATOS_ADD_MODELERS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_MODELERS_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos {

namespace Python {

void  AddCustomModelersToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_MODELERS_TO_PYTHON_H_INCLUDED  defined
