//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"


namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;

    // typedef Properties::Pointer                    PropertiesPointer;
    // typedef Mesh<Node<3>, Properties, Element, Condition>   MeshType;

    // typedef ConstitutiveLaw                  ConstitutiveLawBaseType;
    // typedef ConstitutiveLaw::Pointer          ConstitutiveLawPointer;
    // typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;

    // void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
    // 				      ConstitutiveLawPointer ThisConstitutiveLaw )
    // {
    //   ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
    // }

    void  AddCustomConstitutiveLawsToPython()
    {
      // class_< MaterialsContainer >( "MaterialsContainer", init<>() )
      // 	.def( "PushBack", Push_Back_Constitutive_Laws )
      // 	;
    }
    
  }  // namespace Python.
}  // namespace Kratos.
