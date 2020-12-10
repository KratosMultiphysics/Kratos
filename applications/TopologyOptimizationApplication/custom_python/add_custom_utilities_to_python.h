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

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

  void  AddCustomUtilitiesToPython(pybind11::module& m);

}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
