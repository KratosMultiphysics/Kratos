//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
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
//yield criteria
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
//flow rules
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/bingham_viscoplastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/viscoplastic_flow_rule.hpp"

#include "custom_constitutive/hyperelastic_viscoplastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_viscoplastic_2D_plain_strain_law.hpp"

#include "custom_constitutive/hencky_mc_3D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_2D_law.hpp"

#include "custom_constitutive/hencky_mc_UP_3D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_UP_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

typedef FlowRule::Pointer                        FlowRulePointer;
typedef MPMFlowRule::Pointer                  MPMFlowRulePointer;
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

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    // Hyperelastic Viscoplastic
    class_< HyperElasticViscoplastic3DLaw, typename HyperElasticViscoplastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HyperElasticViscoplastic3DLaw")
    .def(init<>() )
    .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_< HyperElasticViscoplasticPlaneStrain2DLaw, typename HyperElasticViscoplasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HyperElasticViscoplasticPlaneStrain2DLaw")
    .def(init<>() )
    .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    // Hencky Mohr Coulomb
    class_< HenckyMCPlasticPlaneStrain2DLaw, typename HenckyMCPlasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlasticPlaneStrain2DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_< HenckyMCPlastic3DLaw, typename HenckyMCPlastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlastic3DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_< HenckyMCPlasticUP3DLaw, typename HenckyMCPlasticUP3DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlasticUP3DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;

    class_< HenckyMCPlasticPlaneStrainUP2DLaw, typename HenckyMCPlasticPlaneStrainUP2DLaw::Pointer, ConstitutiveLaw >
    (m, "HenckyMCPlasticPlaneStrainUP2DLaw")
    .def(init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
}

}  // namespace Python.
}  // namespace Kratos.
