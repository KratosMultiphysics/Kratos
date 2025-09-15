//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_ADD_CUSTOM_BOUNDING_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_BOUNDING_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

void  AddCustomBoundingToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_BOUNDING_TO_PYTHON_H_INCLUDED  defined
