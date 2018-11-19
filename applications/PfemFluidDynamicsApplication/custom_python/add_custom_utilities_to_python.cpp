//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes

// External includes

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"

// Project includes
#include "includes/node.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

#include "custom_utilities/two_step_v_p_settings.h"

namespace Kratos
{

  namespace Python
  {

    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {

      namespace py = pybind11;
    }

  }  // namespace Python.

} // Namespace Kratos

