//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

#if !defined (KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/flow_rules/non_linear_rate_dependent_plastic_flow_rule.hpp"
#include "custom_constitutive/yield_criteria/mises_huber_thermal_yield_criterion.hpp"
#include "custom_constitutive/hardening_laws/johnson_cook_thermal_hardening_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic-plastic thermal isotropic constitutive law J2 in plane strain 2D
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements
 */

class JohnsonCookThermalPlastic3DLaw : public HyperElasticPlastic3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef ParticleFlowRule::Pointer                FlowRulePointer;
    typedef ParticleYieldCriterion::Pointer    YieldCriterionPointer;
    typedef ParticleHardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of JohnsonCookThermalPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalPlastic3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    JohnsonCookThermalPlastic3DLaw();


    JohnsonCookThermalPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    JohnsonCookThermalPlastic3DLaw(const JohnsonCookThermalPlastic3DLaw& rOther);


    /**
     * Assignment operator.
     */
    JohnsonCookThermalPlastic3DLaw& operator=(const JohnsonCookThermalPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~JohnsonCookThermalPlastic3DLaw() override;

    /**
     * Operators
     */



    /**
     * Calculates the Temperature of the domain (element)
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rTemperature the calculated temperature to be returned
     */
    double& CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
					double & rTemperature) override;

    /**
     * Operations needed by the base class:
     */


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    //int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

    bool Has( const Variable<double>& rThisVariable ) override;

    double & GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

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

    Matrix mElasticLeftCauchyGreen;

    FlowRulePointer       mpFlowRule;

    YieldCriterionPointer mpYieldCriterion;

    HardeningLawPointer   mpHardeningLaw;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElasticPlastic3DLaw);

        rSerializer.save("mElasticLeftCauchyGreen", mElasticLeftCauchyGreen);
        rSerializer.save("mpFlowRule", mpFlowRule);
        rSerializer.save("mpYieldCriterion", mpYieldCriterion);
        rSerializer.save("mpHardeningLaw", mpHardeningLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElasticPlastic3DLaw);

        rSerializer.load("mElasticLeftCauchyGreen", mElasticLeftCauchyGreen);
        rSerializer.load("mpFlowRule", mpFlowRule);
        rSerializer.load("mpYieldCriterion", mpYieldCriterion);
        rSerializer.load("mpHardeningLaw", mpHardeningLaw);
    }



}; // Class JohnsonCookThermalPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED defined
