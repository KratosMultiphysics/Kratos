//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined(KRATOS_POROMECHANICS_ADD_TRILINOS_SCHEMES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_POROMECHANICS_ADD_TRILINOS_SCHEMES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

namespace Kratos {
namespace Python {

void  AddTrilinosSchemesToPython(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_POROMECHANICS_ADD_TRILINOS_SCHEMES_TO_PYTHON_H_INCLUDED 
