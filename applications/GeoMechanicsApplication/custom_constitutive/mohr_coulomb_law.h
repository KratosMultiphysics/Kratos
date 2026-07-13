// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi,
//                   Wijtze Pieter Kikstra
//

#pragma once

#include "custom_constitutive/principal_stresses.hpp"
#include "includes/constitutive_law.h"
#include <memory>

namespace Kratos
{

class ConstitutiveLawDimension;
class CoulombImpl;

class KRATOS_API(GEO_MECHANICS_APPLICATION) MohrCoulombLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombLaw);

    MohrCoulombLaw();
    ~MohrCoulombLaw() override;
    explicit MohrCoulombLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    // Copying is not allowed. Use member `Clone` instead.
    MohrCoulombLaw(const MohrCoulombLaw&)            = delete;
    MohrCoulombLaw& operator=(const MohrCoulombLaw&) = delete;

    // Moving is supported
    MohrCoulombLaw(MohrCoulombLaw&&) noexcept;
    MohrCoulombLaw& operator=(MohrCoulombLaw&&) noexcept;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;
    SizeType                               WorkingSpaceDimension() override;
    bool                                   IsIncremental() override;
    bool                                   RequiresInitializeMaterialResponse() override;
    StressMeasure                          GetStressMeasure() override;
    [[nodiscard]] SizeType                 GetStrainSize() const override;
    StrainMeasure                          GetStrainMeasure() override;
    void    InitializeMaterial(const Properties&     rMaterialProperties,
                               const Geometry<Node>& rElementGeometry,
                               const Vector&         rShapeFunctionsValues) override;
    void    InitializeMaterialResponseCauchy(Parameters& rValues) override;
    void    GetLawFeatures(Features& rFeatures) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    int&    GetValue(const Variable<int>& rVariable, int& rValue) override;
    using ConstitutiveLaw::GetValue;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using ConstitutiveLaw::SetValue;
    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters) override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mStrainVectorFinalized;
    std::unique_ptr<CoulombImpl>              mpCoulombImpl;
    bool                                      mIsModelInitialized = false;
    // Each sub-step's elastic predictor is allowed to overshoot the yield surface by at most
    // this fraction of the stress magnitude
    double      mMaxRelativeOvershoot = 1.0;
    std::size_t mMaxNumberOfSubSteps  = 1;

    [[nodiscard]] Vector CalculateTrialStressVector(const Vector& rStrainVector, const Properties& rProperties) const;

    /// <summary>
    /// Estimates how many strain sub-steps are needed to integrate the constitutive law
    /// accurately. It subdivides the strain increment such that each sub-step overshoots by at most a
    /// target fraction of the stress magnitude. The result is clamped to [1, MaxNumberOfSubSteps].
    /// </summary>
    /// <param name="rImpl"></param>
    /// <param name="rTrialPrincipalStresses"></param>
    /// <param name="rElasticMatrix"></param>
    /// <returns></returns>
    std::size_t CalculateAdaptiveNumberOfSubSteps(std::unique_ptr<CoulombImpl>& rImpl,
                                                  const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                  const Matrix& rElasticMatrix);

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class MohrCoulombLaw

} // namespace Kratos
