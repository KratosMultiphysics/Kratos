//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
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

//hardening laws
#include "custom_constitutive/custom_hardening_laws/hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"

//yield criteria
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"

//flow rules
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"


//constitutive laws
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_constitutive/linear_elastic_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_orthotropic_3D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_U_P_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_U_P_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_plastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_plane_stress_2D_law.hpp"

#include "custom_constitutive/isotropic_damage_simo_ju_3D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef FlowRule::Pointer                        FlowRulePointer;
typedef YieldCriterion::Pointer            YieldCriterionPointer;
typedef HardeningLaw::Pointer                HardeningLawPointer;
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


    //Hyperelastic laws

    class_< HyperElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElastic3DLaw",
      init<>() )
    ;
    
    class_< HyperElasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlaneStrain2DLaw",
      init<>() )
    ;

    class_< HyperElasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticAxisym2DLaw",
      init<>() )
    ;


    //Hyperelastic laws U-P

    class_< HyperElasticUP3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticUP3DLaw",
      init<>() )
    ;


    class_< HyperElasticUPPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticUPPlaneStrain2DLaw",
      init<>() )
    ;

    class_< HyperElasticUPAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticUPAxisym2DLaw",
      init<>() )
    ;


    //Linear Elastic laws

    class_< LinearElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElastic3DLaw",
      init<>() )
    ;

    class_< LinearElasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticPlaneStrain2DLaw",
      init<>() )
    ;

    class_< LinearElasticPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticPlaneStress2DLaw",
      init<>() )
    ;

    class_< LinearElasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticAxisym2DLaw",
      init<>() )
    ;

    class_< LinearElasticOrthotropic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticOrthotropic3DLaw",
      init<>() )
    ;

    //Hyperelastic Plastic laws

    class_<HyperElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlastic3DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HyperElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HyperElasticPlasticAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;


    //Hyperelastic Plastic laws U-P

    class_<HyperElasticPlasticUP3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticUP3DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HyperElasticPlasticUPPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticUPPlaneStrain2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_<HyperElasticPlasticUPAxisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticUPAxisym2DLaw",
      init<>() )
      .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    //Hyperelastic Plastic J2 specilization laws 

    class_<HyperElasticPlasticJ23DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticJ23DLaw",
      init<>() )
    ;

    class_<HyperElasticPlasticJ2PlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticJ2PlaneStrain2DLaw",
      init<>() )
    ;

    class_<HyperElasticPlasticJ2Axisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticJ2Axisym2DLaw",
      init<>() )
    ;

    //Hyperelastic Plastic J2 specilization laws U-P

    class_<HyperElasticPlasticUPJ23DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticUPJ23DLaw",
      init<>() )
    ;

    class_<HyperElasticPlasticUPJ2PlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticUPJ2PlaneStrain2DLaw",
      init<>() )
    ;

    class_<HyperElasticPlasticUPJ2Axisym2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticPlasticUPJ2Axisym2DLaw",
      init<>() )
    ;

    //Linear Elastic Plastic laws 

    class_< LinearElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticPlastic3DLaw",
      init<>() )
    ;

    class_< LinearElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticPlasticPlaneStrain2DLaw",
      init<>() )
    ;

    class_< LinearElasticPlasticPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "LinearElasticPlasticPlaneStress2DLaw",
      init<>() )
    ;

    //Isotropic Damage laws 

    class_< IsotropicDamageSimoJu3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamageSimoJu3DLaw",
      init<>() )
    ;

    class_< IsotropicDamageSimoJuPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamageSimoJuPlaneStrain2DLaw",
      init<>() )
    ;

    class_< IsotropicDamageSimoJuPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "IsotropicDamageSimoJuPlaneStress2DLaw",
      init<>() )
    ;

}

}  // namespace Python.
}  // namespace Kratos.
