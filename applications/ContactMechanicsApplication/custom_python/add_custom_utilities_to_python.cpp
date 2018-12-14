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

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/rigid_body_element_creation_utility.hpp"
#include "custom_utilities/contact_domain_utilities.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
  namespace py = pybind11;

  py::class_<RigidBodyElementCreationUtility>(m,"RigidBodyCreationUtility")
      .def(py::init<>())
      .def("CreateRigidBodyElement",&RigidBodyElementCreationUtility::CreateRigidBodyElement)
      ;
}

}  // namespace Python.

} // Namespace Kratos
