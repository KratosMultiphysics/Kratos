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

#include "custom_constitutive/coulomb_with_tension_cut_off_impl.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

class ConstitutiveLawDimension;

class KRATOS_API(GEO_MECHANICS_APPLICATION) MohrCoulombWithTensionCutOff : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombWithTensionCutOff);

    MohrCoulombWithTensionCutOff() = default;
    explicit MohrCoulombWithTensionCutOff(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    // Copying is not allowed. Use member `Clone` instead.
    MohrCoulombWithTensionCutOff(const MohrCoulombWithTensionCutOff&)            = delete;
    MohrCoulombWithTensionCutOff& operator=(const MohrCoulombWithTensionCutOff&) = delete;

    // Moving is supported
    MohrCoulombWithTensionCutOff(MohrCoulombWithTensionCutOff&&) noexcept            = default;
    MohrCoulombWithTensionCutOff& operator=(MohrCoulombWithTensionCutOff&&) noexcept = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;
    SizeType                               WorkingSpaceDimension() override;
    bool                                   IsIncremental() override;
    bool                                   RequiresInitializeMaterialResponse() override;
    StressMeasure                          GetStressMeasure() override;
    [[nodiscard]] SizeType                 GetStrainSize() const override;
    StrainMeasure                          GetStrainMeasure() override;
    void                                   InitializeMaterial(const Properties&     rMaterialProperties,
                                                              const Geometry<Node>& rElementGeometry,
                                                              const Vector&         rShapeFunctionsValues) override;
    void    InitializeMaterialResponseCauchy(Parameters& rValues) override;
    void    GetLawFeatures(Features& rFeatures) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
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
    CoulombWithTensionCutOffImpl              mCoulombWithTensionCutOffImpl;
    bool                                      mIsModelInitialized = false;

    [[nodiscard]] Vector CalculateTrialStressVector(const Vector& rStrainVector,
                                                    double        YoungsModulus,
                                                    double        PoissonsRatio) const;
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class MohrCoulombWithTensionCutOff

} // namespace Kratos
