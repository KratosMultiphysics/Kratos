//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/continuum_laws/local_damage_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalLocalDamage3DLaw : public LocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    ThermalLocalDamage3DLaw();

    /// Second Constructor
    ThermalLocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /// Copy Constructor
    ThermalLocalDamage3DLaw (const ThermalLocalDamage3DLaw& rOther);

    /// Destructor
    virtual ~ThermalLocalDamage3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    ConstitutiveLaw::Pointer Clone() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Material response entry points.
    //
    // This family is an infinitesimal-strain formulation: for element-provided
    // infinitesimal strain the PK2, Kirchhoff and Cauchy entry points execute
    // the SAME thermo-damage constitutive calculation (total strain minus
    // thermal strain -> elastic predictor -> local-damage return mapping ->
    // damaged stress/tangent). No finite-deformation stress/tangent
    // transformation is performed and no detF scaling is applied.

    void CalculateMaterialResponsePK2 (Parameters & rValues) override;

    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    void FinalizeMaterialResponsePK2 (Parameters & rValues) override;

    void FinalizeMaterialResponseKirchhoff (Parameters & rValues) override;

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Lifecycle
    //
    // The damage state is initialized through InitializeMaterial (flow-rule
    // history/threshold) and subsequently managed through the finalization
    // (commit/restore with IS_CONVERGED). No local-damage state is initialized
    // through the InitializeMaterialResponse hooks, so they are not required.

    bool RequiresInitializeMaterialResponse() override;

    /**
     * @brief Re-establishes the transient material Properties on the flow-rule
     * hardening law from the given current Properties. The hardening-law
     * Properties reference is deliberately NOT serialized (serializer.h:
     * HardeningLaw::save), so after a restart it is NULL until re-established.
     * This method is automatically invoked by every entry point that enters the
     * damage hierarchy (material response, finalization, parameter-aware
     * CalculateValue), so restarted evolution works without any manual repair.
     * Does NOT modify the committed damage/history state.
     */
    void ReinitializeMaterialProperties(const Properties& rMaterialProperties);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Computes the specialized thermo-mechanical vector outputs from the current
     * constitutive state:
     *   THERMAL_STRAIN_VECTOR    = epsilon_th
     *   THERMAL_STRESS_VECTOR    = damage_factor * C * epsilon_th
     *   MECHANICAL_STRESS_VECTOR = damage_factor * C * epsilon_total
     * where damage_factor is the SAME current (1 - d) factor that governs the
     * current total response, so that
     *   stress == MECHANICAL_STRESS_VECTOR - THERMAL_STRESS_VECTOR.
     * The evaluation is read-only: it reads the committed maximum equivalent
     * strain and never commits damage/history or modifies the supplied strain.
     */
    Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * Computes the specialized thermo-mechanical tensor outputs
     * (THERMAL_STRAIN_TENSOR, THERMAL_STRESS_TENSOR, MECHANICAL_STRESS_TENSOR)
     * as the tensor representations of the corresponding vector outputs,
     * obtained by reusing the vector CalculateValue as the single source of
     * truth.
     */
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& CalculateNodalReferenceTemperature ( const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature);

    virtual void CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature);

    /**
     * @brief Validates only the constitutive parameters genuinely consumed by
     * the thermo-damage response/finalization.
     * @note Shape-function derivatives, the deformation gradient and detF are
     * NOT consumed by this infinitesimal-strain implementation, so they are not
     * required (the generic CheckAllParameters requires them unconditionally,
     * which prevented execution under StructuralMechanicsApplication elements).
     */
    void CheckThermalDamageParameters(Parameters& rValues) const;

    /**
     * @brief Common small-strain thermal-damage response shared by the PK2,
     * Kirchhoff and Cauchy entry points:
     *   total infinitesimal strain
     *   -> interpolated temperature / reference temperature
     *   -> thermal strain
     *   -> mechanical strain = total strain - thermal strain (local copy)
     *   -> elastic predictor
     *   -> existing local-damage return mapping
     *   -> existing damage tangent
     *   -> total damaged stress.
     */
    virtual void CalculateThermalDamageResponse(Parameters& rValues);

    /**
     * @brief Common thermal-damage finalization shared by the PK2, Kirchhoff
     * and Cauchy entry points. IS_CONVERGED selects commit (true) or restore
     * (false) of the equilibrium damage/history state.
     */
    virtual void FinalizeThermalDamageResponse(Parameters& rValues);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LocalDamage3DLaw )
        // The stateful damage/history state lives in the flow rule and must be
        // preserved on serialization (restart).
        rSerializer.save("mpFlowRule", mpFlowRule);
        rSerializer.save("mpYieldCriterion", mpYieldCriterion);
        rSerializer.save("mpHardeningLaw", mpHardeningLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LocalDamage3DLaw )
        rSerializer.load("mpFlowRule", mpFlowRule);
        rSerializer.load("mpYieldCriterion", mpYieldCriterion);
        rSerializer.load("mpHardeningLaw", mpHardeningLaw);
    }

}; // Class ThermalLocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined
