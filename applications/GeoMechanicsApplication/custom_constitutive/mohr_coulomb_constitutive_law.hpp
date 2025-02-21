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
#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"
#include "custom_constitutive/coulomb_yield_function.hpp"
#include "includes/serializer.h"

namespace Kratos
{

class ConstitutiveLawDimension;

class KRATOS_API(GEO_MECHANICS_APPLICATION) MohrCoulombConstitutiveLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombConstitutiveLaw);

    // MohrCoulombConstitutiveLaw();
    explicit MohrCoulombConstitutiveLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    ~MohrCoulombConstitutiveLaw() override;

    // Copying is not allowed. Use member `Clone` instead.
    MohrCoulombConstitutiveLaw(const MohrCoulombConstitutiveLaw&)            = delete;
    MohrCoulombConstitutiveLaw& operator=(const MohrCoulombConstitutiveLaw&) = delete;

    // Moving is supported
    MohrCoulombConstitutiveLaw(MohrCoulombConstitutiveLaw&&) noexcept            = default;
    MohrCoulombConstitutiveLaw& operator=(MohrCoulombConstitutiveLaw&&) noexcept = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    using ConstitutiveLaw::GetValue;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    using ConstitutiveLaw::SetValue;
    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMohrCoulomb(const Properties& rProp, Vector& rCautchyStressVector);

    void CalculatePK2Stress(const Vector& rStrainVector, Vector& rStressVector, ConstitutiveLaw::Parameters& rValues);
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues);
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    Vector NormalizeVector(Vector& rVector);
    Matrix ConvertVectorToDiagonalMatrix(Vector& rVector);
    Matrix CalculateRotationMatrix(Matrix& eigenVectorsMatrix);
    void CheckRotationMatrix(Matrix& rRotationMatrix);
    Vector RotatePrincipalStresses(Matrix& rPrincipalStressMatrix, Matrix& rRotationMatrix);
    Vector ReturnStressAtElasticZone(Vector& rTrailStressVector);
    Vector ReturnStressAtAxialZone(Vector& rPrincipalTrialStressVector, double tensionCutoff, Matrix& rRotationMatrix);
    Vector ReturnStressAtCornerReturnZone(Vector& rPrincipalTrialStressVector,
                                          Matrix& rRotationMatrix,
                                          Vector& rCornerPoint);
    Vector ReturnStressAtRegularFailureZone(Vector&                     rPrincipalTrialStressVector,
                                            const CoulombYieldFunction& rCoulombYieldFunction,
                                            Matrix&                     rRotationMatrix,
                                            double                      friction_angle,
                                            double                      cohesion);


    Vector CalculateCornerPoint(double rFrictionAngle, double rCohesion, double rTensionCutoff, double rApex);

    // Member Variables
    double mStateVariable;

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mDeltaStrainVector;
    Vector                                    mStrainVectorFinalized;
    bool                                      mIsCutoffActive = false;

    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

}; // Class MohrCoulombConstitutiveLaw

} // namespace Kratos