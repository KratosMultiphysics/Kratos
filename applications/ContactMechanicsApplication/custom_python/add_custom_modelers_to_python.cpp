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
#include "custom_python/add_custom_modelers_to_python.h"

// Meshers
#include "custom_modelers/contact_domain_2D_modeler.hpp"
#include "custom_modelers/contact_domain_3D_modeler.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomModelersToPython(pybind11::module& m)
{

  using namespace pybind11;

  // Class that allows 3D contact domain spatial search
  class_<ContactDomain3DModeler, typename ContactDomain3DModeler::Pointer, MeshModeler>
      (m,"ContactDomain3DModeler")
      .def(init< >())
      ;

  // Class that allows 2D contact domain spatial search
  class_<ContactDomain2DModeler, typename ContactDomain2DModeler::Pointer, MeshModeler>
      (m,"ContactDomain2DModeler")
      .def(init< >())
      ;
     
}

}  // namespace Python.

} // Namespace Kratos

