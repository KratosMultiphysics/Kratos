//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined (KRATOS_HENCKY_BORJA_CAM_CLAY_PLASTIC_3D_LAW_H_INCLUDED)
#define       KRATOS_HENCKY_BORJA_CAM_CLAY_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hencky_plastic_3d_law.hpp"
#include "custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.hpp"
#include "custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.hpp"
#include "custom_constitutive/hardening_laws/cam_clay_hardening_law.hpp"


namespace Kratos
{
/**
 * Defines a critical-state constitutive law for associative two-invariant Modified Cam Clay model formulated by (Borja, 1998)
 * With stress split in an volumetric and deviatoric parts
 * This material law is defined by the parameters needed by the yield criterion:
 * PRE_CONSOLIDATION_STRESS, OVER_CONSOLIDATION_RATIO, SWELLING_SLOPE, NORMAL_COMPRESSION_SLOPE, CRITICAL_STATE_LINE
 * INITIAL_SHEAR_MODULUS, ALPHA_SHEAR
 * The functionality is designed for large displacements and considering linear ln(v) - ln(p_c) Cam Clay hardening
 * For reference, please refer to: (Borja, 1998)
*/


class HenckyBorjaCamClayPlastic3DLaw
    : public HenckyElasticPlastic3DLaw

{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef MPMFlowRule::Pointer                FlowRulePointer;
    typedef MPMYieldCriterion::Pointer    YieldCriterionPointer;
    typedef MPMHardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of HyperElasticPlasticJ2PlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HenckyBorjaCamClayPlastic3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HenckyBorjaCamClayPlastic3DLaw();


    HenckyBorjaCamClayPlastic3DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    HenckyBorjaCamClayPlastic3DLaw (const HenckyBorjaCamClayPlastic3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticPlasticJ2PlaneStrain2DLaw& operator=(const HyperElasticPlasticJ2PlaneStrain2DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~HenckyBorjaCamClayPlastic3DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo) override;



    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

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

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


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


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }



}; // Class HyperElasticPlasticMohrCoulombPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_HENCKY_BORJA_CAM_CLAY_PLASTIC_3D_LAW_H_INCLUDED defined
