//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Philipp Bucher
//

#if !defined(KRATOS_CONVERGENCE_ACCELERATORS_PYTHON_H_INCLUDED )
#define  KRATOS_CONVERGENCE_ACCELERATORS_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes

namespace Kratos {
namespace Python {

void AddTrilinosConvergenceAcceleratorsToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_CONVERGENCE_ACCELERATORS_PYTHON_H_INCLUDED  defined
