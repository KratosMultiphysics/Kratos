//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/continuum_laws/nonlocal_damage_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalNonlocalDamage3DLaw : public NonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalNonlocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalNonlocalDamage3DLaw();

    /// Second Constructor
    ThermalNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /// Copy Constructor
    ThermalNonlocalDamage3DLaw (const ThermalNonlocalDamage3DLaw& rOther);

    /// Destructor
    virtual ~ThermalNonlocalDamage3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Material response entry points.
    //
    // This family is an infinitesimal-strain formulation: for element-provided
    // infinitesimal strain the PK2, Kirchhoff and Cauchy entry points execute
    // the SAME thermo-NONLOCAL-damage constitutive calculation (total strain
    // minus thermal strain -> elastic predictor -> nonlocal return mapping
    // driven by mNonlocalEquivalentStrain -> damaged stress/tangent). No
    // finite-deformation transformation is performed.

    void CalculateMaterialResponsePK2 (Parameters & rValues) override;

    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    void FinalizeMaterialResponsePK2 (Parameters & rValues) override;

    void FinalizeMaterialResponseKirchhoff (Parameters & rValues) override;

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Lifecycle
    //
    // The damage/history state is initialized through InitializeMaterial and
    // managed through the (stateful) finalization. No local-damage state is
    // initialized through the InitializeMaterialResponse hooks.

    bool RequiresInitializeMaterialResponse() override;

    /**
     * @brief Re-establishes the transient material Properties on the flow-rule
     * hardening law after a serialization/restart (the transient Properties are
     * not serialized). Required before a restored law can produce a further
     * material response. Does NOT modify the committed damage/history state.
     */
    void ReinitializeMaterialProperties(const Properties& rMaterialProperties);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Parameter-aware scalar outputs. For LOCAL_EQUIVALENT_STRAIN this
    // recomputes the local driving quantity from the CURRENT kinematics
    // supplied by the element and stores it in the flow-rule state exposed by
    // GetValue(LOCAL_EQUIVALENT_STRAIN). Has(LOCAL_EQUIVALENT_STRAIN) is NOT
    // implemented, so the generic StructuralMechanics
    // CalculateOnConstitutiveLaw path reaches this method.

    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& CalculateNodalReferenceTemperature ( const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature);

    virtual void CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature);

    /**
     * @brief Validates only the parameters genuinely consumed by the thermal
     * nonlocal formulation (no ShapeFunctionsDerivatives / F / detF unless a
     * specific path uses them).
     */
    void CheckThermalNonlocalDamageParameters(Parameters& rValues) const;

    /**
     * @brief Common LOCAL equivalent-strain calculation:
     *   total strain -> thermal strain -> mechanical strain (local copy)
     *   -> elastic predictor -> existing CalculateLocalReturnMapping,
     * updating the flow-rule thermal variable exposed by
     * GetValue(LOCAL_EQUIVALENT_STRAIN). It never commits the damage/history
     * state.
     * @param rStressVector Optional output (legacy flag path requests stress).
     * @param rConstitutiveMatrix Optional output (legacy flag path requests tangent).
     */
    void CalculateLocalEquivalentStrain(
        Parameters& rValues,
        Vector* pStressVector = nullptr,
        Matrix* pConstitutiveMatrix = nullptr);

    /**
     * @brief Common small-strain thermal NONLOCAL response. When
     * INITIALIZE_MATERIAL_RESPONSE is present (temporary legacy compatibility)
     * only the LOCAL driving-quantity calculation is executed.
     */
    virtual void CalculateThermalNonlocalDamageResponse(Parameters& rValues);

    /**
     * @brief Common thermal NONLOCAL finalization (IS_CONVERGED commit/restore),
     * driven by mNonlocalEquivalentStrain.
     */
    virtual void FinalizeThermalNonlocalDamageResponse(Parameters& rValues);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonlocalDamage3DLaw )
        // Stateful nonlocal damage state must survive serialization (restart).
        rSerializer.save("mpFlowRule", mpFlowRule);
        rSerializer.save("mpYieldCriterion", mpYieldCriterion);
        rSerializer.save("mpHardeningLaw", mpHardeningLaw);
        rSerializer.save("mNonlocalEquivalentStrain", mNonlocalEquivalentStrain);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonlocalDamage3DLaw )
        rSerializer.load("mpFlowRule", mpFlowRule);
        rSerializer.load("mpYieldCriterion", mpYieldCriterion);
        rSerializer.load("mpHardeningLaw", mpHardeningLaw);
        rSerializer.load("mNonlocalEquivalentStrain", mNonlocalEquivalentStrain);
    }

}; // Class ThermalNonlocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined
