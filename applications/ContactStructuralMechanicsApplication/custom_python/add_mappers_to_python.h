// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 


#if !defined(KRATOS_ADD_MAPPERS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_MAPPERS_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

void  AddCustomMappersToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_MAPPERS_TO_PYTHON_H_INCLUDED  defined 
