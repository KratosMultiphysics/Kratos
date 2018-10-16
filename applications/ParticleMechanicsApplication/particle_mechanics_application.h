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


#if !defined(KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "solid_mechanics_application.h"
#include "particle_mechanics_application_variables.h"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/kratos_application.h"

#include "includes/condition.h"  
#include "includes/ublas_interface.h"

#include "containers/flags.h"

/* CONDITIONS */
#include "custom_conditions/mpm_base_load_condition.h"
#include "custom_conditions/mpm_point_load_condition.h"
#include "custom_conditions/mpm_line_load_condition_2d.h"
#include "custom_conditions/mpm_surface_load_condition_3d.h"

//---element
#include "custom_elements/updated_lagrangian.hpp"
#include "custom_elements/updated_lagrangian_UP.hpp"
#include "custom_elements/updated_lagrangian_quadrilateral.hpp"

//---constitutive laws
#include "custom_constitutive/hyperelastic_viscoplastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_viscoplastic_2D_plain_strain_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_UP_2D_law.hpp"
#include "custom_constitutive/hencky_mc_3D_law.hpp"
#include "custom_constitutive/hencky_mc_UP_3D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_3D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_3D_law.hpp"

//---flow rules
#include "custom_constitutive/flow_rules/viscoplastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/bingham_viscoplastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/mc_plastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/mc_strain_softening_plastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.hpp"

//---yield criteria
#include "custom_constitutive/yield_criteria/mc_yield_criterion.hpp"
#include "custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.hpp"

//---hardening laws
#include "custom_constitutive/hardening_laws/exponential_strain_softening_law.hpp"
#include "custom_constitutive/hardening_laws/cam_clay_hardening_law.hpp"

namespace Kratos
{

/// Short class definition.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for particle mechanics problems. 
 * Currently developed methods are: (1) Material Point Method
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

    // Elements
    const UpdatedLagrangian mUpdatedLagrangian2D3N;
    const UpdatedLagrangian mUpdatedLagrangian3D4N;
    const UpdatedLagrangianUP mUpdatedLagrangianUP2D3N;
    //const UpdatedLagrangianUP mUpdatedLagrangianUP3D4N;
    const UpdatedLagrangianQuadrilateral mUpdatedLagrangian2D4N;
    const UpdatedLagrangianQuadrilateral mUpdatedLagrangian3D8N;
    //const UpdatedLagrangianUPQuadrilateral mUpdatedLagrangianUP2D4N;
    //const TotalLagrangian mTotalLagrangian2D3N;
    //const TotalLagrangian mTotalLagrangian3D4N;
    
    // Conditions
    const MPMPointLoadCondition mMPMPointLoadCondition2D1N;
    const MPMPointLoadCondition mMPMPointLoadCondition3D1N;
    const MPMLineLoadCondition2D mMPMLineLoadCondition2D2N;
    const MPMSurfaceLoadCondition3D mMPMSurfaceLoadCondition3D3N;
    const MPMSurfaceLoadCondition3D mMPMSurfaceLoadCondition3D4N;

    // Constitutive laws
    // CL: Hyperelastic ViscoPlastic laws
    const HyperElasticViscoplastic3DLaw                     mHyperElasticViscoplastic3DLaw;
    const HyperElasticViscoplasticPlaneStrain2DLaw          mHyperElasticViscoplasticPlaneStrain2DLaw;
    // CL: Mohr Coulomb
    const HenckyMCPlastic3DLaw                              mHenckyMCPlastic3DLaw;
    const HenckyMCPlasticPlaneStrain2DLaw                   mHenckyMCPlasticPlaneStrain2DLaw;
    const HenckyMCPlasticUP3DLaw                            mHenckyMCPlasticUP3DLaw;
    const HenckyMCPlasticPlaneStrainUP2DLaw                 mHenckyMCPlasticPlaneStrainUP2DLaw;
    // CL: Mohr Coulomb Strain Softening
    const HenckyMCStrainSofteningPlastic3DLaw               mHenckyMCStrainSofteningPlastic3DLaw;
    const HenckyMCStrainSofteningPlasticPlaneStrain2DLaw    mHenckyMCStrainSofteningPlasticPlaneStrain2DLaw;
    // CL: Borja Cam Clay
    const HenckyBorjaCamClayPlastic3DLaw                    mHenckyBorjaCamClayPlastic3DLaw;
    const HenckyBorjaCamClayPlasticPlaneStrain2DLaw         mHenckyBorjaCamClayPlasticPlaneStrain2DLaw;

    // Flow Rules
    const ViscoplasticFlowRule                      mViscoplasticFlowRule;
    const BinghamViscoplasticFlowRule               mBinghamViscoplasticFlowRule;
    const MCPlasticFlowRule                         mMCPlasticFlowRule;
    const MCStrainSofteningPlasticFlowRule          mMCStrainSofteningPlasticFlowRule;
    const BorjaCamClayPlasticFlowRule               mBorjaCamClayPlasticFlowRule;
    //const NonLinearAssociativePlasticFlowRule     mNonLinearAssociativePlasticFlowRule;
    //const LinearAssociativePlasticFlowRule        mLinearAssociativePlasticFlowRule;
    //const IsotropicDamageFlowRule                 mIsotropicDamageFlowRule;
    //const DruckerPragerFlowRule                   mDruckerPragerFlowRule;

    // Yield Criteria
    const MCYieldCriterion                          mMCYieldCriterion;
    const ModifiedCamClayYieldCriterion             mModifiedCamClayYieldCriterion;
    //const MisesHuberYieldCriterion                mMisesHuberYieldCriterion;
    //const SimoJuYieldCriterion                    mSimoJuYieldCriterion;
    //const DruckerPragerYieldCriterion             mDruckerPragerYieldCriterion;

    // Hardening Laws
    const ExponentialStrainSofteningLaw             mExponentialStrainSofteningLaw;
    const CamClayHardeningLaw                       mCamClayHardeningLaw;
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


