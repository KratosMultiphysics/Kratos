//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_FLUID_DYNAMICS_ADD_TRILINOS_PROCESSES_TO_PYTHON_H_INCLUDED
#define KRATOS_FLUID_DYNAMICS_ADD_TRILINOS_PROCESSES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

namespace Kratos {
namespace Python {

void AddTrilinosProcessesToPython(pybind11::module& m);

}
}

#endif // KRATOS_FLUID_DYNAMICS_ADD_TRILINOS_PROCESSES_TO_PYTHON_H_INCLUDED