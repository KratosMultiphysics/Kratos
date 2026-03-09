//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

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

    void AddCustomProcessesToPython(pybind11::module& m);

  } // namespace Python
} // namespace Kratos
