//
//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:            Duan Wenjie $
//   Date:                $Date:                March 2016 $
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

using namespace boost::python;

typedef FlowRule::Pointer                        FlowRulePointer;
typedef MPMFlowRule::Pointer                        MPMFlowRulePointer;
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
    class_<HyperElasticViscoplastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticViscoplastic3DLaw",
      init<>() )
    .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<HyperElasticViscoplasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HyperElasticViscoplasticPlaneStrain2DLaw",
      init<>() )
    .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<HenckyMCPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMCPlasticPlaneStrain2DLaw",
      init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<HenckyMCPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMCPlastic3DLaw",
      init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<HenckyMCPlasticUP3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMCPlasticUP3DLaw",
      init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
    class_<HenckyMCPlasticPlaneStrainUP2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "HenckyMCPlasticPlaneStrainUP2DLaw",
      init<>() )
    .def( init<MPMFlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
    ;
}

}  // namespace Python.
}  // namespace Kratos.
