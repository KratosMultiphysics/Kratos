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

class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceCoulombWithTensionCutOff : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCoulombWithTensionCutOff);

    InterfaceCoulombWithTensionCutOff() = default;
    explicit InterfaceCoulombWithTensionCutOff(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    // Copying is not allowed. Use member `Clone` instead.
    InterfaceCoulombWithTensionCutOff(const InterfaceCoulombWithTensionCutOff&)            = delete;
    InterfaceCoulombWithTensionCutOff& operator=(const InterfaceCoulombWithTensionCutOff&) = delete;

    // Moving is supported
    InterfaceCoulombWithTensionCutOff(InterfaceCoulombWithTensionCutOff&&) noexcept = default;
    InterfaceCoulombWithTensionCutOff& operator=(InterfaceCoulombWithTensionCutOff&&) noexcept = default;

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
    void    InitializeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    using ConstitutiveLaw::GetValue;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using ConstitutiveLaw::SetValue;
    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;
    void    CalculateMaterialResponseCauchy(Parameters& rConstitutiveLawParameters) override;
    void    FinalizeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters) override;
    Matrix& CalculateValue(Parameters&             rConstitutiveLawParameters,
                           const Variable<Matrix>& rVariable,
                           Matrix&                 rValue) override;
    using ConstitutiveLaw::CalculateValue;

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mTractionVector;
    Vector                                    mTractionVectorFinalized;
    Vector                                    mRelativeDisplacementVectorFinalized;
    CoulombWithTensionCutOffImpl              mCoulombWithTensionCutOffImpl;
    bool                                      mIsModelInitialized = false;

    [[nodiscard]] Geo::SigmaTau CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                             double NormalStiffness,
                                                             double ShearStiffness) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class InterfaceMohrCoulombWithTensionCutOff

} // namespace Kratos
