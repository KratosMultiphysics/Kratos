//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"

#include "custom_utilities/two_step_v_p_settings.h"

namespace Kratos
{

  namespace Python
  {

    void  AddCustomUtilitiesToPython(pybind11::module& m);

  }  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined 
