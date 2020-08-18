//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#ifdef KRATOS_PYTHON
#include <pybind11/pybind11.h>

#include "add_trilinos_strategies_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosParticleMechanicsTrilinosExtension,m) {
    AddTrilinosStrategiesToPython(m);
}

}// namespace Python
}// namespace Kratos

#endif