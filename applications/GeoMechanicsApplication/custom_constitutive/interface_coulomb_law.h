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

#include "includes/constitutive_law.h"
#include "sigma_tau.hpp"
#include <memory>

namespace Kratos
{

class ConstitutiveLawDimension;
class CoulombImpl;

class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceCoulombLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCoulombLaw);

    InterfaceCoulombLaw();
    ~InterfaceCoulombLaw() override;
    explicit InterfaceCoulombLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    // Copying is not allowed. Use member `Clone` instead.
    InterfaceCoulombLaw(const InterfaceCoulombLaw&)            = delete;
    InterfaceCoulombLaw& operator=(const InterfaceCoulombLaw&) = delete;

    // Moving is supported
    InterfaceCoulombLaw(InterfaceCoulombLaw&&) noexcept;
    InterfaceCoulombLaw& operator=(InterfaceCoulombLaw&&) noexcept;

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
    void    InitializeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    int&    GetValue(const Variable<int>& rVariable, int& rValue) override;
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
    std::unique_ptr<CoulombImpl>              mpCoulombImpl;
    bool                                      mIsModelInitialized = false;

    [[nodiscard]] Geo::SigmaTau CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                             double NormalStiffness,
                                                             double ShearStiffness) const;
    std::size_t CalculateAdaptiveNumberOfSubSteps(const Geo::SigmaTau& rTrialTraction, const Matrix& rElasticMatrix);

    double      mMaxRelativeOvershoot       = 1.0;
    int         mMaxNumberOfSubSteps        = 1;
    std::size_t mCalculatedNumberOfSubSteps = 0;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class InterfaceMohrCoulomb

} // namespace Kratos
