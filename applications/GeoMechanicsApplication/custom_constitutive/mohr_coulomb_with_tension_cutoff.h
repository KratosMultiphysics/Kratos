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

// System includes
#include <cmath>

// Project includes
#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/tension_cutoff.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
class ConstitutiveLawDimension;
class Serializer;

class KRATOS_API(GEO_MECHANICS_APPLICATION) MohrCoulombWithTensionCutOff : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombWithTensionCutOff);

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
    SizeType                               GetStrainSize() const override;
    StrainMeasure                          GetStrainMeasure() override;
    void InitializeMaterial(const Properties& rMaterialProperties,
                                    const Geometry<Node>& rElementGeometry,
                                    const Vector& rShapeFunctionsValues) override;
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    void SetValue(const Variable<Vector >& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    int  Check(const Properties&   rMaterialProperties,
               const GeometryType& rElementGeometry,
               const ProcessInfo&  rCurrentProcessInfo) const override;
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters) override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    Vector RotatePrincipalStresses(const Vector& rPrincipalStressVector, const Matrix& rRotationMatrix) const;

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mStrainVectorFinalized;
    CoulombYieldSurface                       mCoulombYieldSurface;
    TensionCutoff                             mTensionCutOff;

    void CheckProperty(const Properties& rMaterialProperties, const Kratos::Variable<double>& rVariable) const;
    void   CalculateTrialStressVector(const Vector&                rStrainVector,
                                      Vector&                      rStressVector,
                                      double YoungModulus,
                                                              double PoissonRatio) const;
    double CalculateApex(double FrictionAngle, double Cohesion) const;
    Vector CalculateCornerPoint(double FrictionAngle, double Cohesion, double TensionCutoff) const;
    Vector ReturnStressAtAxialZone(const Vector& rPrincipalTrialStressVector, double TensionCutoff) const;
    Vector ReturnStressAtCornerReturnZone(const Vector& rPrincipalTrialStressVector,
                                          const Vector& rCornerPoint) const;
    Vector ReturnStressAtRegularFailureZone(const Vector&              rPrincipalTrialStressVector,
                                            const CoulombYieldSurface& rCoulombYieldFunction,
                                            double                     FrictionAngle,
                                            double                     Cohesion) const;

    // Serialization
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class MohrCoulombWithTensionCutOff
} // namespace Kratos
