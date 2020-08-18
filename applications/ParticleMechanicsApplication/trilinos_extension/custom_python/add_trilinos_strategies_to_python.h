//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Manuel Messmer
//

#ifndef KRATOS_ADD_TRILINOS_STRATEGIES_TO_PYTHON_H_INCLUDED
#define KRATOS_ADD_TRILINOS_STRATEGIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

namespace Kratos {
namespace Python {

void AddTrilinosStrategiesToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // KRATOS_ADD_TRILINOS_STRATEGIES_TO_PYTHON_H_INCLUDED