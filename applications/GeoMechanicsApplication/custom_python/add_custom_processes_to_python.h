// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_ADD_CUSTOM_PROCESSES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_CUSTOM_PROCESSES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>
// Project includes
#include "includes/define_python.h"

namespace Kratos::Python
{
void AddCustomProcessesToPython(pybind11::module& m);
} // namespace Kratos::Python

#endif // KRATOS_ADD_CUSTOM_PROCESSES_TO_PYTHON_H_INCLUDED  defined
