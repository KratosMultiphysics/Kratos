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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "processes/process.h"
#include "custom_python/add_custom_modelers_to_python.h"

// Meshers
#include "custom_modelers/contact_domain_2D_modeler.hpp"
//#include "custom_modelers/contact_domain_3D_modeler.hpp"

// Bounding Boxes


namespace Kratos
{

namespace Python
{

  typedef MeshModeler                        MeshModelerBaseType;

  void  AddCustomModelersToPython()
  {

    using namespace boost::python;
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)

    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)


    //class that allows 2D contact domain spatial search
    class_<ContactDomain2DModeler, bases<MeshModelerBaseType>, boost::noncopyable >
      ("ContactDomain2DModeler", init< >())
      ;
     
  }

}  // namespace Python.

} // Namespace Kratos

