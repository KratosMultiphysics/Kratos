// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//


#if !defined(KRATOS_GEO_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_GEO_MECHANICS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes
#include "geo_mechanics_application_variables.h"

// conditions
#include "custom_conditions/U_Pw_force_condition.hpp"
#include "custom_conditions/U_Pw_face_load_condition.hpp"
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_FIC_condition.hpp"
#include "custom_conditions/U_Pw_face_load_interface_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_interface_condition.hpp"
#include "custom_conditions/line_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp"
#include "custom_conditions/surface_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.hpp"

// elements
#include "custom_elements/transient_Pw_element.hpp"
#include "custom_elements/steady_state_Pw_element.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"
#include "custom_elements/U_Pw_small_strain_link_interface_element.hpp"
#include "custom_elements/U_Pw_small_strain_FIC_element.hpp"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_elements/drained_U_Pw_small_strain_element.hpp"
#include "custom_elements/undrained_U_Pw_small_strain_element.hpp"
#include "custom_elements/U_Pw_updated_lagrangian_element.hpp"
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.hpp"
#include "custom_elements/U_Pw_updated_lagrangian_FIC_element.hpp"
#include "custom_elements/small_strain_U_Pw_diff_order_axisymmetric_element.hpp"
#include "custom_elements/U_Pw_small_strain_axisymmetric_element.hpp"
#include "custom_elements/U_Pw_small_strain_axisymmetric_FIC_element.hpp"

/* geo structural element */
#include "custom_elements/geo_cr_beam_element_3D2N.hpp"
#include "custom_elements/geo_cr_beam_element_2D2N.hpp"
#include "custom_elements/geo_cr_beam_element_linear_2D2N.hpp"
#include "custom_elements/geo_cr_beam_element_linear_3D2N.hpp"
#include "custom_elements/geo_truss_element.hpp"
#include "custom_elements/geo_linear_truss_element.hpp"
#include "custom_elements/geo_cable_element.hpp"

// constitutive models
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"
#include "custom_constitutive/elastic_isotropic_K0_3d_law.h"
#include "custom_constitutive/linear_elastic_plane_strain_K0_law.h"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.h"

#include "custom_constitutive/small_strain_udsm_3D_law.hpp"
#include "custom_constitutive/small_strain_udsm_2D_plane_strain_law.hpp"
#include "custom_constitutive/small_strain_udsm_2D_interface_law.hpp"
#include "custom_constitutive/small_strain_udsm_3D_interface_law.hpp"

#include "custom_constitutive/small_strain_umat_3D_law.hpp"
#include "custom_constitutive/small_strain_umat_2D_plane_strain_law.hpp"
#include "custom_constitutive/small_strain_umat_2D_interface_law.hpp"
#include "custom_constitutive/small_strain_umat_3D_interface_law.hpp"

#include "custom_constitutive/linear_elastic_2D_interface_law.h"
#include "custom_constitutive/linear_elastic_3D_interface_law.h"


namespace Kratos {

///@name Kratos Globals
///@{

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
class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoMechanicsApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosGeoMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosGeoMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosGeoMechanicsApplication();

    /// Destructor.
    virtual ~KratosGeoMechanicsApplication(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register() override;



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
    virtual std::string Info() const override
    {
        return "KratosGeoMechanicsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
          KRATOS_WATCH("in KratosGeoMechanicsApplication");
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

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    // const Elem2D   mElem2D;
    // const Elem3D   mElem3D;

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

    // elements
    // transient one-phase flow elements:
    const TransientPwElement<2,3> mTransientPwElement2D3N;
    const TransientPwElement<2,4> mTransientPwElement2D4N;
    const TransientPwElement<3,4> mTransientPwElement3D4N;
    const TransientPwElement<3,8> mTransientPwElement3D8N;

    const TransientPwElement<2,6> mTransientPwElement2D6N;
    const TransientPwElement<2,8> mTransientPwElement2D8N;
    const TransientPwElement<2,9> mTransientPwElement2D9N;
    const TransientPwElement<3,10> mTransientPwElement3D10N;
    const TransientPwElement<3,20> mTransientPwElement3D20N;
    const TransientPwElement<3,27> mTransientPwElement3D27N;

    // Steady-State one-phase flow elements:
    const SteadyStatePwElement<2,3> mSteadyStatePwElement2D3N;
    const SteadyStatePwElement<2,4> mSteadyStatePwElement2D4N;
    const SteadyStatePwElement<3,4> mSteadyStatePwElement3D4N;
    const SteadyStatePwElement<3,8> mSteadyStatePwElement3D8N;

    const SteadyStatePwElement<2,6> mSteadyStatePwElement2D6N;
    const SteadyStatePwElement<2,8> mSteadyStatePwElement2D8N;
    const SteadyStatePwElement<2,9> mSteadyStatePwElement2D9N;
    const SteadyStatePwElement<3,10> mSteadyStatePwElement3D10N;
    const SteadyStatePwElement<3,20> mSteadyStatePwElement3D20N;
    const SteadyStatePwElement<3,27> mSteadyStatePwElement3D27N;

    // small strain elements:
    const UPwSmallStrainElement<2,3> mUPwSmallStrainElement2D3N;
    const UPwSmallStrainElement<2,4> mUPwSmallStrainElement2D4N;
    const UPwSmallStrainElement<3,4> mUPwSmallStrainElement3D4N;
    const UPwSmallStrainElement<3,8> mUPwSmallStrainElement3D8N;

    const UPwSmallStrainElement<2,6> mUPwSmallStrainElement2D6N;
    const UPwSmallStrainElement<2,8> mUPwSmallStrainElement2D8N;
    const UPwSmallStrainElement<2,9> mUPwSmallStrainElement2D9N;
    const UPwSmallStrainElement<3,10> mUPwSmallStrainElement3D10N;
    const UPwSmallStrainElement<3,20> mUPwSmallStrainElement3D20N;
    const UPwSmallStrainElement<3,27> mUPwSmallStrainElement3D27N;

    // small strain drained elements:
    const DrainedUPwSmallStrainElement<2,3> mDrainedUPwSmallStrainElement2D3N;
    const DrainedUPwSmallStrainElement<2,4> mDrainedUPwSmallStrainElement2D4N;
    const DrainedUPwSmallStrainElement<3,4> mDrainedUPwSmallStrainElement3D4N;
    const DrainedUPwSmallStrainElement<3,8> mDrainedUPwSmallStrainElement3D8N;

    // small strain undrained elements:
    const UndrainedUPwSmallStrainElement<2,3> mUndrainedUPwSmallStrainElement2D3N;
    const UndrainedUPwSmallStrainElement<2,4> mUndrainedUPwSmallStrainElement2D4N;
    const UndrainedUPwSmallStrainElement<3,4> mUndrainedUPwSmallStrainElement3D4N;
    const UndrainedUPwSmallStrainElement<3,8> mUndrainedUPwSmallStrainElement3D8N;

    // FIC elements
    const UPwSmallStrainFICElement<2,3> mUPwSmallStrainFICElement2D3N;
    const UPwSmallStrainFICElement<2,4> mUPwSmallStrainFICElement2D4N;
    const UPwSmallStrainFICElement<3,4> mUPwSmallStrainFICElement3D4N;
    const UPwSmallStrainFICElement<3,8> mUPwSmallStrainFICElement3D8N;

    // Small strain different order elements
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D6N;
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D8N;
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D9N;
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D10N;
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D20N;
    const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D27N;

    // small strain axisymmtric elements:
    const UPwSmallStrainAxisymmetricElement<2,3> mUPwSmallStrainAxisymmetricElement2D3N;
    const UPwSmallStrainAxisymmetricElement<2,4> mUPwSmallStrainAxisymmetricElement2D4N;
    const UPwSmallStrainAxisymmetricElement<2,6> mUPwSmallStrainAxisymmetricElement2D6N;
    const UPwSmallStrainAxisymmetricElement<2,8> mUPwSmallStrainAxisymmetricElement2D8N;
    const UPwSmallStrainAxisymmetricElement<2,9> mUPwSmallStrainAxisymmetricElement2D9N;

    const UPwSmallStrainAxisymmetricFICElement<2,3> mUPwSmallStrainAxisymmetricFICElement2D3N;
    const UPwSmallStrainAxisymmetricFICElement<2,4> mUPwSmallStrainAxisymmetricFICElement2D4N;

    const SmallStrainUPwDiffOrderAxisymmetricElement mSmallStrainUPwDiffOrderAxisymmetricElement2D6N;
    const SmallStrainUPwDiffOrderAxisymmetricElement mSmallStrainUPwDiffOrderAxisymmetricElement2D8N;
    const SmallStrainUPwDiffOrderAxisymmetricElement mSmallStrainUPwDiffOrderAxisymmetricElement2D9N;

    // interface elements
    const UPwSmallStrainInterfaceElement<2,4> mUPwSmallStrainInterfaceElement2D4N;
    const UPwSmallStrainInterfaceElement<3,6> mUPwSmallStrainInterfaceElement3D6N;
    const UPwSmallStrainInterfaceElement<3,8> mUPwSmallStrainInterfaceElement3D8N;

    const UPwSmallStrainLinkInterfaceElement<2,4> mUPwSmallStrainLinkInterfaceElement2D4N;
    const UPwSmallStrainLinkInterfaceElement<3,6> mUPwSmallStrainLinkInterfaceElement3D6N;
    const UPwSmallStrainLinkInterfaceElement<3,8> mUPwSmallStrainLinkInterfaceElement3D8N;

    // Updated-Lagrangian elements:
    const UPwUpdatedLagrangianElement<2,3> mUPwUpdatedLagrangianElement2D3N;
    const UPwUpdatedLagrangianElement<2,4> mUPwUpdatedLagrangianElement2D4N;
    const UPwUpdatedLagrangianElement<3,4> mUPwUpdatedLagrangianElement3D4N;
    const UPwUpdatedLagrangianElement<3,8> mUPwUpdatedLagrangianElement3D8N;

    const UPwUpdatedLagrangianElement<2,6> mUPwUpdatedLagrangianElement2D6N;
    const UPwUpdatedLagrangianElement<2,8> mUPwUpdatedLagrangianElement2D8N;
    const UPwUpdatedLagrangianElement<2,9> mUPwUpdatedLagrangianElement2D9N;
    const UPwUpdatedLagrangianElement<3,10> mUPwUpdatedLagrangianElement3D10N;
    const UPwUpdatedLagrangianElement<3,20> mUPwUpdatedLagrangianElement3D20N;
    const UPwUpdatedLagrangianElement<3,27> mUPwUpdatedLagrangianElement3D27N;

    const UPwUpdatedLagrangianFICElement<2,3> mUPwUpdatedLagrangianFICElement2D3N;
    const UPwUpdatedLagrangianFICElement<2,4> mUPwUpdatedLagrangianFICElement2D4N;
    const UPwUpdatedLagrangianFICElement<3,4> mUPwUpdatedLagrangianFICElement3D4N;
    const UPwUpdatedLagrangianFICElement<3,8> mUPwUpdatedLagrangianFICElement3D8N;

    // Updated-Lagrangian different order elements
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D6N;
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D8N;
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement2D9N;
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement3D10N;
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement3D20N;
    const UpdatedLagrangianUPwDiffOrderElement mUpdatedLagrangianUPwDiffOrderElement3D27N;

    // geo structural element
    const GeoCrBeamElement2D2N mGeoCrBeamElement2D2N;
    const GeoCrBeamElement3D2N mGeoCrBeamElement3D2N;
    const GeoCrBeamElementLinear2D2N mGeoCrBeamElementLinear2D2N;
    const GeoCrBeamElementLinear3D2N mGeoCrBeamElementLinear3D2N;
    const GeoTrussElement<2,2> mGeoTrussElement2D2N;
    const GeoTrussElement<3,2> mGeoTrussElement3D2N;

    const GeoLinearTrussElement<2,2> mGeoLinearTrussElement2D2N;
    const GeoLinearTrussElement<3,2> mGeoLinearTrussElement3D2N;

    const GeoCableElement<2,2> mGeoCableElement2D2N;
    const GeoCableElement<3,2> mGeoCableElement3D2N;

    // conditions
    const UPwForceCondition<2,1> mUPwForceCondition2D1N;
    const UPwForceCondition<3,1> mUPwForceCondition3D1N;
    const UPwFaceLoadCondition<2,2> mUPwFaceLoadCondition2D2N;
    const UPwFaceLoadCondition<3,3> mUPwFaceLoadCondition3D3N;
    const UPwFaceLoadCondition<3,4> mUPwFaceLoadCondition3D4N;
    const UPwNormalFaceLoadCondition<2,2> mUPwNormalFaceLoadCondition2D2N;
    const UPwNormalFaceLoadCondition<3,3> mUPwNormalFaceLoadCondition3D3N;
    const UPwNormalFaceLoadCondition<3,4> mUPwNormalFaceLoadCondition3D4N;
    const UPwNormalFluxCondition<2,2> mUPwNormalFluxCondition2D2N;
    const UPwNormalFluxCondition<3,3> mUPwNormalFluxCondition3D3N;
    const UPwNormalFluxCondition<3,4> mUPwNormalFluxCondition3D4N;

    const UPwFaceLoadCondition<2,3> mUPwFaceLoadCondition2D3N;

    const UPwFaceLoadInterfaceCondition<2,2> mUPwFaceLoadInterfaceCondition2D2N;
    const UPwFaceLoadInterfaceCondition<3,4> mUPwFaceLoadInterfaceCondition3D4N;
    const UPwNormalFluxInterfaceCondition<2,2> mUPwNormalFluxInterfaceCondition2D2N;
    const UPwNormalFluxInterfaceCondition<3,4> mUPwNormalFluxInterfaceCondition3D4N;

    const UPwNormalFluxFICCondition<2,2> mUPwNormalFluxFICCondition2D2N;
    const UPwNormalFluxFICCondition<3,3> mUPwNormalFluxFICCondition3D3N;
    const UPwNormalFluxFICCondition<3,4> mUPwNormalFluxFICCondition3D4N;

    const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D3N;
    const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D3N;
    const LineNormalFluidFlux2DDiffOrderCondition mLineNormalFluidFluxDiffOrderCondition2D3N;
    const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D6N;
    const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D8N;
    const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D9N;
    const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D6N;
    const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D8N;
    const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D9N;
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D6N;
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D8N;
    const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D9N;

    // constitutive models
    const BilinearCohesive3DLaw             mBilinearCohesive3DLaw;
    const BilinearCohesive2DLaw             mBilinearCohesive2DLaw;
    const LinearPlaneStrainK0Law            mLinearPlaneStrainK0Law;
    const GeoLinearElasticPlaneStrain2DLaw  mLinearElasticPlaneStrain2DLaw;
    const ElasticIsotropicK03DLaw           mElasticIsotropicK03DLaw;
    const GeoLinearElasticPlaneStress2DLaw  mLinearElasticPlaneStress2DLaw;

    const SmallStrainUDSM3DLaw            mSmallStrainUDSM3DLaw;
    const SmallStrainUDSM2DPlaneStrainLaw mSmallStrainUDSM2DPlaneStrainLaw;
    const SmallStrainUDSM2DInterfaceLaw   mSmallStrainUDSM2DInterfaceLaw;
    const SmallStrainUDSM3DInterfaceLaw   mSmallStrainUDSM3DInterfaceLaw;

    const SmallStrainUMAT3DLaw            mSmallStrainUMAT3DLaw;
    const SmallStrainUMAT2DPlaneStrainLaw mSmallStrainUMAT2DPlaneStrainLaw;
    const SmallStrainUMAT2DInterfaceLaw   mSmallStrainUMAT2DInterfaceLaw;
    const SmallStrainUMAT3DInterfaceLaw   mSmallStrainUMAT3DInterfaceLaw;

    const LinearElastic2DInterfaceLaw     mLinearElastic2DInterfaceLaw;
    const LinearElastic3DInterfaceLaw     mLinearElastic3DInterfaceLaw;

    /// Assignment operator.
    KratosGeoMechanicsApplication& operator=(KratosGeoMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosGeoMechanicsApplication(KratosGeoMechanicsApplication const& rOther);


    ///@}

}; // Class KratosGeoMechanicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_GEO_MECHANICS_APPLICATION_H_INCLUDED  defined
