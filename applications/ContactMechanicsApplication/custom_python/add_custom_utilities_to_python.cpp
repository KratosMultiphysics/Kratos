//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Utilities
#include "custom_utilities/rigid_body_element_creation_utility.hpp"
#include "custom_utilities/contact_domain_utilities.hpp"


namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;
	
	class_< RigidBodyElementCreationUtility, boost::noncopyable > 
	  ("RigidBodyCreationUtility", init<>())
	  .def("CreateRigidBodyElement",&RigidBodyElementCreationUtility::CreateRigidBodyElement)
	  ;

  }



}  // namespace Python.

} // Namespace Kratos
