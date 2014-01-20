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


//yield criteria

//flow rule

//constitutive laws
//#include "custom_constitutive/hencky_plastic_3d_law.hpp"
//#include "custom_constitutive/hencky_plastic_plane_strain_2d_law.hpp"
//#include "custom_constitutive/hencky_plastic_axisym_2d_law.hpp"
#include "custom_constitutive/hencky_cam_clay_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_hencky_cam_clay_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_cam_clay_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_J2_axisym_2D_law.hpp"

namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;

	typedef FlowRule::Pointer                        FlowRulePointer;
	typedef YieldCriterion::Pointer            YieldCriterionPointer;
	typedef HardeningLaw::Pointer                HardeningLawPointer;
	typedef ConstitutiveLaw                  ConstitutiveLawBaseType;

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

    class_<NonLinearHenckyCamClayPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "NonLinearHenckyCamClayPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<LinearHenckyCamClayPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearHenckyCamClayPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<NonLinearHenckyCamClayPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "NonLinearHenckyCamClayPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyJ2PlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyJ2PlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
       
    class_<HenckyJ2PlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyJ2PlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
      // 	;
    }
    
  }  // namespace Python.
}  // namespace Kratos.
