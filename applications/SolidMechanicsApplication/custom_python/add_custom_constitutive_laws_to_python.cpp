//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
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
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_2D_law.hpp"
#include "custom_constitutive/linear_elastic_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"



namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;

    typedef Properties::Pointer                    PropertiesPointer;
    typedef Mesh<Node<3>, Properties, Element, Condition>   MeshType;

    typedef ConstitutiveLaw                  ConstitutiveLawBaseType;
    typedef ConstitutiveLaw::Pointer          ConstitutiveLawPointer;
    typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;

    void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
				      ConstitutiveLawPointer ThisConstitutiveLaw )
    {
      ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
    }

    void  AddCustomConstitutiveLawsToPython()
    {
      class_< MaterialsContainer >( "MaterialsContainer", init<>() )
	.def( "PushBack", Push_Back_Constitutive_Laws )
	;


      class_< HyperElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
	( "HyperElastic3DLaw",
	  init<>() )
	;


     class_< LinearElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
	( "LinearElastic3DLaw",
	  init<>() )
	;


     class_< LinearElastic2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
	( "LinearElastic2DLaw",
	  init<>() )
	;

     class_< LinearElasticPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
	( "LinearElasticPlaneStress2DLaw",
	  init<>() )
	;

      class_< HyperElastic2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
	( "HyperElastic2DLaw",
	  init<>() )
	;

    }
    
  }  // namespace Python.
}  // namespace Kratos.
