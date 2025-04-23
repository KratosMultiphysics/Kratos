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
#include <optional>

// Project includes
#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/tension_cutoff.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

class ConstitutiveLawDimension;

class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceMohrCoulombWithTensionCutOff : public ConstitutiveLaw
{
public:
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
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    using ConstitutiveLaw::GetValue;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using ConstitutiveLaw::SetValue;
    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters) override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

private:
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mStrainVectorFinalized;
    CoulombYieldSurface                       mCoulombYieldSurface;
    TensionCutoff                             mTensionCutOff;
    bool                                      mIsModelInitialized = false;

    [[nodiscard]] Vector CalculateTrialStressVector(const Vector& rStrainVector,
                                                    double        YoungsModulus,
                                                    double        PoissonsRatio) const;
    [[nodiscard]] Matrix MakeConstitutiveMatrix(double NormalStiffness, double ShearStiffness) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class InterfaceMohrCoulombWithTensionCutOff

} // namespace Kratos
