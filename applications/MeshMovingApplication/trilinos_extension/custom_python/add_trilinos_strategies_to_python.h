//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_MESH_MOVING_ADD_TRILINOS_STRATEGIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_MESH_MOVING_ADD_TRILINOS_STRATEGIES_TO_PYTHON_H_INCLUDED

// System includes


// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"


namespace Kratos {
namespace Python {

void  AddMeshMovingStrategies(pybind11::module& m);

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_MESH_MOVING_ADD_TRILINOS_STRATEGIES_TO_PYTHON_H_INCLUDED defined
