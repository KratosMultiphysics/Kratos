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


#if !defined(KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include "solid_mechanics_application.h"
//#include "structural_application.h"
//#include "pfem_solid_mechanics_application.h"
// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"

#include "includes/condition.h"   //#include "costum_conditions/nuova_condizione.h nel caso non sia gia stata definita
#include "includes/ublas_interface.h"

#include "containers/flags.h"

//condition
//#include "custom_conditions/mpm_line_load_2D_condition.hpp"
//#include "custom_conditions/mpm_line_load_3D_condition.hpp"


//element
#include "custom_elements/updated_lagrangian.hpp"
#include "custom_elements/updated_lagrangian_UP.hpp"
#include "custom_elements/updated_lagrangian_quadrilateral.hpp"
//#include "custom_elements/updated_lagrangian_UP_quadrilateral.hpp"
//#include "custom_elements/total_lagrangian.hpp"

//constitutive laws
#include "custom_constitutive/hyperelastic_viscoplastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_viscoplastic_2D_plain_strain_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_UP_2D_law.hpp"


#include "custom_constitutive/hencky_mc_3D_law.hpp"
#include "custom_constitutive/hencky_mc_UP_3D_law.hpp"
//flow rules
#include "custom_constitutive/flow_rules/viscoplastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/bingham_viscoplastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"
//#include "custom_constitutive/flow_rules/drucker_prager_flow_rule.hpp"
//yield criteria
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
//#include "custom_constitutive/yield_criteria/drucker_prager_yield_criterion.hpp"
namespace Kratos
{
///@name Type Definitions
///@{
typedef array_1d<double,3> Vector3;
typedef array_1d<double,6> Vector6;
///@}

///@name Kratos Globals
///@{

// Variables definition

//element
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, COUNTER )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MP_NUMBER )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MP_BOOL )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MP_MATERIAL_ID )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, WEIGHT )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_MASS )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DENSITY )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_VOLUME )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_KINETIC_ENERGY )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_STRAIN_ENERGY )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_TOTAL_ENERGY )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_PRESSURE )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_JACOBIAN )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_EQUIVALENT_PLASTIC_STRAIN )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_CONSTITUTIVE_PRESSURE )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, NODAL_MPRESSURE )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, AUX_PRESSURE)
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, AUX_MP_PRESSURE)


//constitutive law

KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )

KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION,  double, DILATANCY_COEFFICIENT )


KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, COHESION )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_DILATANCY_ANGLE )

////nodal dofs
//KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, IMPOSED_DISPLACEMENT )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, AUX_R)
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, AUX_T)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_R_VEL )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_T_VEL )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_R_ACC )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_T_ACC )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, NODAL_LUMPED_MASS)

KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_VELOCITY )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_ACCELERATION )

//MP element variable
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, GAUSS_COORD )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_DISPLACEMENT )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_VELOCITY )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_ACCELERATION )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_MP_VELOCITY )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, AUX_MP_ACCELERATION )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_VOLUME_ACCELERATION )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, MP_CAUCHY_STRESS_VECTOR )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, MP_ALMANSI_STRAIN_VECTOR )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, PREVIOUS_MP_CAUCHY_STRESS_VECTOR )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, PREVIOUS_MP_ALMANSI_STRAIN_VECTOR )
KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Matrix, MP_CONSTITUTIVE_MATRIX )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, DISPLACEMENT_AUX )
//grid node variable
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_MOMENTUM )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_INERTIA )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_INTERNAL_FORCE )


///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KratosParticleMechanicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosParticleMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosParticleMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosParticleMechanicsApplication();

    /// Destructor.
    ~KratosParticleMechanicsApplication() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosParticleMechanicsApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in my application");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{



    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{
    const UpdatedLagrangian mUpdatedLagrangian2D3N;
    const UpdatedLagrangian mUpdatedLagrangian3D4N;
    const UpdatedLagrangianUP mUpdatedLagrangianUP2D3N;
    //const UpdatedLagrangianUP mUpdatedLagrangianUP3D4N;
    const UpdatedLagrangianQuadrilateral mUpdatedLagrangian2D4N;
    //const UpdatedLagrangianUPQuadrilateral mUpdatedLagrangianUP2D4N;
    //const TotalLagrangian mTotalLagrangian2D3N;
    //const TotalLagrangian mTotalLagrangian3D4N;
    
    //conditions
    //const MPMLineLoad2DCondition                mMPMLineLoadCondition2D2N;
    //const MPMLineLoad2DCondition                mMPMLineLoadCondition2D3N;
    //const MPMLineLoad3DCondition                mMPMLineLoadCondition3D2N;
    //const MPMLineLoad3DCondition                mMPMLineLoadCondition3D3N;

    //constitutive laws
    const HyperElasticViscoplastic3DLaw                mHyperElasticViscoplastic3DLaw;
    const HyperElasticViscoplasticPlaneStrain2DLaw     mHyperElasticViscoplasticPlaneStrain2DLaw;
    const HenckyMCPlastic3DLaw                         mHenckyMCPlastic3DLaw;
    const HenckyMCPlasticPlaneStrain2DLaw              mHenckyMCPlasticPlaneStrain2DLaw;
    const HenckyMCPlasticUP3DLaw                       mHenckyMCPlasticUP3DLaw;
    const HenckyMCPlasticPlaneStrainUP2DLaw            mHenckyMCPlasticPlaneStrainUP2DLaw;

    //Flow Rules
    const ViscoplasticFlowRule                      mViscoplasticFlowRule;
    const BinghamViscoplasticFlowRule               mBinghamViscoplasticFlowRule;
    const MCPlasticFlowRule                         mMCPlasticFlowRule;
    //const NonLinearAssociativePlasticFlowRule     mNonLinearAssociativePlasticFlowRule;
    //const LinearAssociativePlasticFlowRule        mLinearAssociativePlasticFlowRule;
    //const IsotropicDamageFlowRule                 mIsotropicDamageFlowRule;
    //const DruckerPragerFlowRule                   mDruckerPragerFlowRule;

    //Yield Criteria
    const MCYieldCriterion                          mMCYieldCriterion;
    //const MisesHuberYieldCriterion                mMisesHuberYieldCriterion;
    //const SimoJuYieldCriterion                    mSimoJuYieldCriterion;
    //const DruckerPragerYieldCriterion             mDruckerPragerYieldCriterion;

    //Hardening Laws
    //const NonLinearIsotropicKinematicHardeningLaw mNonLinearIsotropicKinematicHardeningLaw;
    //const LinearIsotropicKinematicHardeningLaw    mLinearIsotropicKinematicHardeningLaw;
    //const ExponentialDamageHardeningLaw           mExponentialDamageHardeningLaw;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosParticleMechanicsApplication& operator=(KratosParticleMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosParticleMechanicsApplication(KratosParticleMechanicsApplication const& rOther);


    ///@}

}; // Class KratosParticleMechanicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED  defined 


