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
// //

#pragma once

// System includes
#include <cmath>

// Project includes
#include "custom_constitutive/coulomb_yield_surface.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"

namespace Kratos
{
class ConstitutiveLawDimension;

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

    int    Check(const Properties&   rMaterialProperties,
                 const GeometryType& rElementGeometry,
                 const ProcessInfo&  rCurrentProcessInfo) const override;
    void   CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& parameters) override;
    void   FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    Vector NormalizeVector(const Vector& rVector) const;
    Matrix ConvertVectorToDiagonalMatrix(const Vector& rVector) const;
    Matrix CalculateRotationMatrix(const Matrix& eigenVectorsMatrix) const;
    Vector RotatePrincipalStresses(const Vector& rPrincipalStressVector, const Matrix& rRotationMatrix) const;

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mStrainVectorFinalized;

    void   CalculateTrialStressVector(const Vector&                rStrainVector,
                                      Vector&                      rStressVector,
                                      ConstitutiveLaw::Parameters& rValues) const;
    Matrix CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues) const;
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
}; // Class MohrCoulombConstitutiveLaw
} // namespace Kratos
