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
#include "custom_python/add_custom_meshers_to_python.h"

// Meshers
#include "custom_meshers/contact_domain_2D_mesher.hpp"
#include "custom_meshers/contact_domain_3D_mesher.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomMeshersToPython(pybind11::module& m)
{

  using namespace pybind11;

  // Class that allows 3D contact domain spatial search
  class_<ContactDomain3DMesher, typename ContactDomain3DMesher::Pointer, Mesher>
      (m,"ContactDomain3DMesher")
      .def(init< >())
      ;

  // Class that allows 2D contact domain spatial search
  class_<ContactDomain2DMesher, typename ContactDomain2DMesher::Pointer, Mesher>
      (m,"ContactDomain2DMesher")
      .def(init< >())
      ;

}

}  // namespace Python.

} // Namespace Kratos

