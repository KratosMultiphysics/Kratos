//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes

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

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
//---yield criteria
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
#include "custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.hpp"

//---hardening laws
#include "custom_constitutive/hardening_laws/exponential_strain_softening_law.hpp"
#include "custom_constitutive/hardening_laws/cam_clay_hardening_law.hpp"

//---flow rules
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/mc_strain_softening_plastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.hpp"

//---constitutive laws
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_UP_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_UP_2D_law.hpp"
#include "custom_constitutive/hencky_mc_3D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mc_UP_3D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_UP_2D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_3D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_3D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_plane_strain_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

typedef MPMFlowRule::Pointer                  MPMFlowRulePointer;
typedef MPMYieldCriterion::Pointer      MPMYieldCriterionPointer;
typedef MPMHardeningLaw::Pointer          MPMHardeningLawPointer;
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

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    // Hyperelastic laws
    class_< HyperElastic3DLaw, typename HyperElastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HyperElastic3DLaw").def(init<>() )
    ;

    class_< HyperElasticPlaneStrain2DLaw, typename HyperElasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HyperElasticPlaneStrain2DLaw").def(init<>() )
    ;

    class_< HyperElasticUP3DLaw, typename HyperElasticUP3DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticUP3DLaw").def(init<>() )
      ;

    class_< HyperElasticPlaneStrainUP2DLaw, typename HyperElasticPlaneStrainUP2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlaneStrainUP2DLaw").def(init<>() )
      ;

    // Hencky Mohr Coulomb
    class_< HenckyMCPlasticPlaneStrain2DLaw, typename HenckyMCPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlasticPlaneStrain2DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    class_< HenckyMCPlastic3DLaw, typename HenckyMCPlastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlastic3DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    class_< HenckyMCPlasticUP3DLaw, typename HenckyMCPlasticUP3DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlasticUP3DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    class_< HenckyMCPlasticPlaneStrainUP2DLaw, typename HenckyMCPlasticPlaneStrainUP2DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlasticPlaneStrainUP2DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    // Hencky Mohr Coulomb Strain Softening
    class_< HenckyMCStrainSofteningPlasticPlaneStrain2DLaw, typename HenckyMCStrainSofteningPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCStrainSofteningPlasticPlaneStrain2DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    class_< HenckyMCStrainSofteningPlastic3DLaw, typename HenckyMCStrainSofteningPlastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCStrainSofteningPlastic3DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    // Hencky Borja Cam Clay
    class_< HenckyBorjaCamClayPlasticPlaneStrain2DLaw, typename HenckyBorjaCamClayPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyBorjaCamClayPlasticPlaneStrain2DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;

    class_< HenckyBorjaCamClayPlastic3DLaw, typename HenckyBorjaCamClayPlastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyBorjaCamClayPlastic3DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, MPMYieldCriterionPointer, MPMHardeningLawPointer>() )
    ;
}

}  // namespace Python.
}  // namespace Kratos.
