//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes
#include "includes/define.h"
#include "includes/define_python.h"

// Project includes

namespace Kratos
{
  namespace Python
  {

    void AddCustomUtilitiesToPython(pybind11::module& m);

  } // namespace Python
} // namespace Kratos

#endif // KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED defined
