//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Last modified by:    $Author:               JMCarbonell $
//   Date:                $Date:              September 2015 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes 
#include <boost/python.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes


namespace Kratos
{
	
  namespace Python
  {

  	
    void  AddCustomProcessesToPython()
    {

      using namespace boost::python;
      typedef Process                                         ProcessBaseType;

    }
 
  }  // namespace Python.

} // Namespace Kratos

