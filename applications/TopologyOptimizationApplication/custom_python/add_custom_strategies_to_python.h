// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_STRATEGIES_PYTHON_H_INCLUDED )
#define  KRATOS_STRATEGIES_PYTHON_H_INCLUDED

// System includes 
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

void  AddCustomStrategiesToPython(pybind11::module& m);
/* { */
/*   using namespace boost::python; */
/* } */

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_STRATEGIES_PYTHON_H_INCLUDED  defined
