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


#if !defined(KRATOS_MPM_APPLICATION_H_INCLUDED )
#define  KRATOS_MPM_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "mpm_application_variables.h"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/kratos_application.h"

#include "includes/condition.h"
#include "includes/ublas_interface.h"

#include "containers/flags.h"

/* CONDITIONS */
#include "custom_conditions/grid_based_conditions/mpm_grid_base_load_condition.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_point_load_condition.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_axisym_point_load_condition.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_line_load_condition_2d.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_axisym_line_load_condition_2d.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_surface_load_condition_3d.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_dirichlet_condition.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.h"

//---element
#include "custom_elements/mpm_updated_lagrangian.hpp"
#include "custom_elements/mpm_updated_lagrangian_UP.hpp"
#include "custom_elements/mpm_updated_lagrangian_PQ.hpp"

//---constitutive laws
#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_axisym_2D_law.hpp"
#include "custom_constitutive/johnson_cook_thermal_plastic_3D_law.hpp"
#include "custom_constitutive/johnson_cook_thermal_plastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/johnson_cook_thermal_plastic_axisym_2D_law.hpp"
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_axisym_2D_law.hpp"
#include "custom_constitutive/hyperelastic_UP_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_UP_2D_law.hpp"
#include "custom_constitutive/hencky_mc_3D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mc_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_mc_UP_3D_law.hpp"
#include "custom_constitutive/hencky_mc_plane_strain_UP_2D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_3D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_mc_strain_softening_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_3D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_borja_cam_clay_axisym_2D_law.hpp"
#include "custom_constitutive/displacement_newtonian_fluid_3D_law.hpp"
#include "custom_constitutive/displacement_newtonian_fluid_plane_strain_2D_law.hpp"

//---flow rules
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
 * for MPM problems.
 * Currently developed methods are: (1) Material Point Method
 */
class KRATOS_API(MPM_APPLICATION) KratosMPMApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMPMApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMPMApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMPMApplication();

    /// Destructor.
    ~KratosMPMApplication() override {}


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
        return "KratosMPMApplication";
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
    const MPMUpdatedLagrangian mMPMUpdatedLagrangian;
    const MPMUpdatedLagrangianUP mMPMUpdatedLagrangianUP;
    const MPMUpdatedLagrangianPQ mMPMUpdatedLagrangianPQ;

    // Deprecated Elements
    const MPMUpdatedLagrangian mMPMUpdatedLagrangian2D3N;
    const MPMUpdatedLagrangian mMPMUpdatedLagrangian3D4N;
    const MPMUpdatedLagrangian mMPMUpdatedLagrangianUP2D3N;
    const MPMUpdatedLagrangian mMPMUpdatedLagrangian2D4N;
    const MPMUpdatedLagrangian mMPMUpdatedLagrangian3D8N;
    const MPMUpdatedLagrangian mMPMUpdatedLagrangianAxisymmetry2D3N;
    const MPMUpdatedLagrangian mMPMUpdatedLagrangianAxisymmetry2D4N;

    // Conditions
    // Grid Conditions:
    const MPMGridPointLoadCondition mMPMGridPointLoadCondition2D1N;
    const MPMGridPointLoadCondition mMPMGridPointLoadCondition3D1N;
    const MPMGridAxisymPointLoadCondition mMPMGridAxisymPointLoadCondition2D1N;
    const MPMGridLineLoadCondition2D mMPMGridLineLoadCondition2D2N;
    const MPMGridAxisymLineLoadCondition2D mMPMGridAxisymLineLoadCondition2D2N;
    const MPMGridSurfaceLoadCondition3D mMPMGridSurfaceLoadCondition3D3N;
    const MPMGridSurfaceLoadCondition3D mMPMGridSurfaceLoadCondition3D4N;
    // MPM Conditions:
    const MPMParticlePenaltyDirichletCondition mMPMParticlePenaltyDirichletCondition;
    const MPMParticlePointLoadCondition mMPMParticlePointLoadCondition;

    // Deprecated Conditions
    const MPMParticlePenaltyDirichletCondition mMPMParticlePenaltyDirichletCondition2D3N;
    const MPMParticlePenaltyDirichletCondition mMPMParticlePenaltyDirichletCondition2D4N;
    const MPMParticlePenaltyDirichletCondition mMPMParticlePenaltyDirichletCondition3D4N;
    const MPMParticlePenaltyDirichletCondition mMPMParticlePenaltyDirichletCondition3D8N;
    const MPMParticlePointLoadCondition mMPMParticlePointLoadCondition2D3N;
    const MPMParticlePointLoadCondition mMPMParticlePointLoadCondition3D4N;
    const MPMParticlePointLoadCondition mMPMParticlePointLoadCondition2D4N;
    const MPMParticlePointLoadCondition mMPMParticlePointLoadCondition3D8N;


    // Constitutive laws
    // CL: Linear Elastic laws
    const LinearElastic3DLaw                                mLinearElastic3DLaw;
    const LinearElasticPlaneStress2DLaw                     mLinearElasticPlaneStress2DLaw;
    const LinearElasticPlaneStrain2DLaw                     mLinearElasticPlaneStrain2DLaw;
    const LinearElasticAxisym2DLaw                          mLinearElasticAxisym2DLaw;
    // CL: Johnson Cooker Thermal Plastic laws
    const JohnsonCookThermalPlastic3DLaw                    mJohnsonCookThermalPlastic3DLaw;
    const JohnsonCookThermalPlastic2DPlaneStrainLaw         mJohnsonCookThermalPlastic2DPlaneStrainLaw;
    const JohnsonCookThermalPlastic2DAxisymLaw              mJohnsonCookThermalPlastic2DAxisymLaw;
    // CL: Hyperelastic laws
    const HyperElastic3DLaw                                 mHyperElastic3DLaw;
    const HyperElasticPlaneStrain2DLaw                      mHyperElasticPlaneStrain2DLaw;
    const HyperElasticAxisym2DLaw                           mHyperElasticAxisym2DLaw;
    const HyperElasticUP3DLaw                               mHyperElasticUP3DLaw;
    const HyperElasticPlaneStrainUP2DLaw                    mHyperElasticPlaneStrainUP2DLaw;
    // CL: Mohr Coulomb
    const HenckyMCPlastic3DLaw                              mHenckyMCPlastic3DLaw;
    const HenckyMCPlasticPlaneStrain2DLaw                   mHenckyMCPlasticPlaneStrain2DLaw;
    const HenckyMCPlasticAxisym2DLaw                        mHenckyMCPlasticAxisym2DLaw;
    const HenckyMCPlasticUP3DLaw                            mHenckyMCPlasticUP3DLaw;
    const HenckyMCPlasticPlaneStrainUP2DLaw                 mHenckyMCPlasticPlaneStrainUP2DLaw;
    // CL: Mohr Coulomb Strain Softening
    const HenckyMCStrainSofteningPlastic3DLaw               mHenckyMCStrainSofteningPlastic3DLaw;
    const HenckyMCStrainSofteningPlasticPlaneStrain2DLaw    mHenckyMCStrainSofteningPlasticPlaneStrain2DLaw;
    const HenckyMCStrainSofteningPlasticAxisym2DLaw         mHenckyMCStrainSofteningPlasticAxisym2DLaw;
    // CL: Borja Cam Clay
    const HenckyBorjaCamClayPlastic3DLaw                    mHenckyBorjaCamClayPlastic3DLaw;
    const HenckyBorjaCamClayPlasticPlaneStrain2DLaw         mHenckyBorjaCamClayPlasticPlaneStrain2DLaw;
    const HenckyBorjaCamClayPlasticAxisym2DLaw              mHenckyBorjaCamClayPlasticAxisym2DLaw;
    // CL: Displacement-based Newtonian Fluid
    const DispNewtonianFluid3DLaw                           mDispNewtonianFluid3DLaw;
    const DispNewtonianFluidPlaneStrain2DLaw                mDispNewtonianFluidPlaneStrain2DLaw;
    // Flow Rules
    const MCPlasticFlowRule                         mMCPlasticFlowRule;
    const MCStrainSofteningPlasticFlowRule          mMCStrainSofteningPlasticFlowRule;
    const BorjaCamClayPlasticFlowRule               mBorjaCamClayPlasticFlowRule;

    // Yield Criteria
    const MCYieldCriterion                          mMCYieldCriterion;
    const ModifiedCamClayYieldCriterion             mModifiedCamClayYieldCriterion;

    // Hardening Laws
    const ExponentialStrainSofteningLaw             mExponentialStrainSofteningLaw;
    const CamClayHardeningLaw                       mCamClayHardeningLaw;

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
    KratosMPMApplication& operator=(KratosMPMApplication const& rOther);

    /// Copy constructor.
    KratosMPMApplication(KratosMPMApplication const& rOther);


    ///@}

}; // Class KratosMPMApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_MPM_APPLICATION_H_INCLUDED  defined


