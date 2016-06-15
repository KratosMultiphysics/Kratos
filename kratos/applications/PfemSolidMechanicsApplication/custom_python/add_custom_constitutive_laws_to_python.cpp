//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
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

//hardening laws
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"

//yield criteria
#include "custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/J2_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/tresca_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/mohr_coulomb_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/matsuoka_nakai_yield_criterion.hpp"

//flow rules
#include "custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"
//#include "custom_constitutive/custom_flow_rules/cam_clay_explicit_plastic_flow_rule.hpp"
//#include "custom_constitutive/custom_flow_rules/linear_cam_clay_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/borja_cam_clay_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/mohr_coulomb_explicit_plastic_flow_rule.hpp"

//#include "custom_constitutive/custom_flow_rules/non_associative_plastic_flow_rule.hpp"
//#include "custom_constitutive/custom_flow_rules/matsuoka_nakai_flow_rule.hpp"



//constitutive laws
//#include "custom_constitutive/hencky_plastic_3d_law.hpp"
//#include "custom_constitutive/hencky_plastic_plane_strain_2d_law.hpp"
//#include "custom_constitutive/hencky_plastic_axisym_2d_law.hpp"
//#include "custom_constitutive/hencky_matsuoka_plane_strain_2D_law.hpp"
//#include "custom_constitutive/hencky_matsuoka_axisym_2D_law.hpp"
//#include "custom_constitutive/hencky_von_eekelen_axisym_2D_law.hpp"
//#include "custom_constitutive/hencky_cam_clay_plane_strain_2D_law.hpp"
//#include "custom_constitutive/hencky_cam_clay_axisym_2D_law.hpp"
//#include "custom_constitutive/linear_hencky_cam_clay_plane_strain_2D_law.hpp"
//#include "custom_constitutive/linear_hencky_cam_clay_axisym_2D_law.hpp"

#include "custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.hpp"
#include "custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.hpp"

#include "custom_constitutive/hencky_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mohr_coulomb_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_mohr_coulomb_plane_strain_2D_law.hpp"

#include "custom_constitutive/hencky_U_P_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.hpp"

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
/*    class_<HenckyElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyElasticPlastic3DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;


    class_<HenckyElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyElasticPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
*/
/*    class_<HenckyElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyElasticPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

*/
/*    class_<HenckyMatsuokaPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMatsuokaPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyMatsuokaPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMatsuokaPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyVonEekelenPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyVonEekelenPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
*/

   /* class_<NonLinearHenckyCamClayPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
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
    class_<LinearHenckyCamClayPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearHenckyCamClayPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;*/

    class_<BorjaHenckyCamClayPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "BorjaHenckyCamClayPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<BorjaHenckyCamClayPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "BorjaHenckyCamClayPlasticPlaneStrain2DLaw",
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

    class_<HenckyPlasticUPJ2Axisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyPlasticUPJ2Axisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyPlasticUPJ2PlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyPlasticUPJ2PlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyPlasticUPTrescaAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyPlasticUPTrescaAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyPlasticUPTrescaPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyPlasticUPTrescaPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HenckyTrescaPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyTrescaPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
      // 	;
    class_<HenckyTrescaPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyTrescaPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<HenckyMohrCoulombPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMohrCoulombPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
      // 	;
    class_<HenckyMohrCoulombPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMohrCoulombPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    }
    
  }  // namespace Python.
}  // namespace Kratos.
