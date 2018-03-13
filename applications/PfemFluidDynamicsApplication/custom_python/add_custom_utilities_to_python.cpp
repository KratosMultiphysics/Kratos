//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/two_step_v_p_settings.h"

namespace Kratos
{
	
  namespace Python
  {
    
    void  AddCustomUtilitiesToPython()
    {

      using namespace boost::python;
    }

  }  // namespace Python.

} // Namespace Kratos

