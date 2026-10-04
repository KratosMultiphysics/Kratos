// License: BSD License
// Kratos default license: kratos/license.txt

#pragma once

#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @brief Small-strain, associated J2 plasticity with sigma_y = K*(epsilon_0+p)^n.
 * @details The committed history is plastic strain and accumulated equivalent
 * plastic strain p. Strain vectors use [xx, yy, zz, 2xy, 2yz, 2xz]; stress
 * vectors use [xx, yy, zz, xy, yz, xz]. As in SmallStrainJ2Plasticity3D,
 * Delta epsilon_p = DeltaGamma*N, N = s_trial/||s_trial||, N:N = 1,
 * and Delta p = sqrt(2/3)*DeltaGamma. Thus DeltaGamma is the multiplier
 * for ||s||-sqrt(2/3)*sigma_y; the multiplier for sigma_eq-sigma_y is Delta p.
 * CalculateMaterialResponse only evaluates trial states; Finalize commits.
 * All stress-measure entry points coincide in this small-strain formulation.
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SmallStrainSwiftJ2Plasticity3D
    : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainSwiftJ2Plasticity3D);
    using BaseType = ConstitutiveLaw;
    using BaseType::Has;
    using BaseType::GetValue;
    using BaseType::SetValue;
    using BaseType::CalculateValue;

    SmallStrainSwiftJ2Plasticity3D() = default;
    SmallStrainSwiftJ2Plasticity3D(const SmallStrainSwiftJ2Plasticity3D&) = default;
    ~SmallStrainSwiftJ2Plasticity3D() override = default;

    ConstitutiveLaw::Pointer Clone() const override;
    SizeType WorkingSpaceDimension() override { return 3; }
    SizeType GetStrainSize() const override { return 6; }
    StressMeasure GetStressMeasure() override { return StressMeasure_Cauchy; }
    void GetLawFeatures(Features& rFeatures) override;
    std::string Info() const override { return "SmallStrainSwiftJ2Plasticity3DLaw"; }

    bool RequiresInitializeMaterialResponse() override { return false; }
    bool RequiresFinalizeMaterialResponse() override { return true; }
    void InitializeMaterial(const Properties& rProperties, const GeometryType& rGeometry,
        const Vector& rShapeFunctionsValues) override;

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;
    void CalculateMaterialResponsePK1(Parameters& rValues) override { CalculateMaterialResponseCauchy(rValues); }
    void CalculateMaterialResponsePK2(Parameters& rValues) override { CalculateMaterialResponseCauchy(rValues); }
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override { CalculateMaterialResponseCauchy(rValues); }

    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;
    void FinalizeMaterialResponsePK1(Parameters& rValues) override { FinalizeMaterialResponseCauchy(rValues); }
    void FinalizeMaterialResponsePK2(Parameters& rValues) override { FinalizeMaterialResponseCauchy(rValues); }
    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override { FinalizeMaterialResponseCauchy(rValues); }

    bool Has(const Variable<double>& rVariable) override;
    bool Has(const Variable<Vector>& rVariable) override;
    double& GetValue(const Variable<double>& rVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    void SetValue(const Variable<double>& rVariable, const double& rValue,
        const ProcessInfo& rProcessInfo) override;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue,
        const ProcessInfo& rProcessInfo) override;

    /// YIELD_STRESS is the Swift flow stress at the committed p, not an input parameter.
    double& CalculateValue(Parameters& rValues, const Variable<double>& rVariable,
        double& rValue) override;

    int Check(const Properties& rProperties, const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo) const override;

private:
    Vector mPlasticStrain = ZeroVector(6);
    double mAccumulatedPlasticStrain = 0.0;

    static void CheckMaterialProperties(const Properties& rProperties);
    static void EvaluateHardening(const Properties& rProperties, const double P,
        double& rYieldStress, double& rHardeningModulus);
    void IntegrateResponse(Parameters& rValues, Vector& rPlasticStrain, double& rP);

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
