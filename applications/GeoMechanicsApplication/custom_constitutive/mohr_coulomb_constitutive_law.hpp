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

        explicit MohrCoulombConstitutiveLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

        // Copying is not allowed. Use member `Clone` instead.
        MohrCoulombConstitutiveLaw(const MohrCoulombConstitutiveLaw&) = delete;
        MohrCoulombConstitutiveLaw& operator=(const MohrCoulombConstitutiveLaw&) = delete;

        // Moving is supported
        MohrCoulombConstitutiveLaw(MohrCoulombConstitutiveLaw&&) noexcept = default;
        MohrCoulombConstitutiveLaw& operator=(MohrCoulombConstitutiveLaw&&) noexcept = default;

        [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

        int Check(const Properties& rMaterialProperties,
                  const GeometryType& rElementGeometry,
                  const ProcessInfo& rCurrentProcessInfo) const override;

        using ConstitutiveLaw::GetValue;
        double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

        using ConstitutiveLaw::SetValue;
        void SetValue(const Variable<double>& rVariable, const double& rValue,
                      const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateMohrCoulomb(ConstitutiveLaw::Parameters& parameters);

        void CalculateTrialStressVector(const Vector& rStrainVector, Vector& rStressVector,
                                        ConstitutiveLaw::Parameters& rValues);
        Matrix CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues) const;
        void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

        Vector NormalizeVector(const Vector& rVector) const;
        Matrix ConvertVectorToDiagonalMatrix(const Vector& rVector) const;
        Matrix CalculateRotationMatrix(const Matrix& eigenVectorsMatrix);
        Vector RotatePrincipalStresses(const Vector& rPrincipalStressVector, Matrix& rRotationMatrix);
        Vector ReturnStressAtElasticZone(const Vector& rTrailStressVector) const;
        Vector ReturnStressAtAxialZone(const Vector& rPrincipalTrialStressVector, const double TensionCutoff) const;
        Vector ReturnStressAtCornerReturnZone(const Vector& rPrincipalTrialStressVector,
                                              const Vector& rCornerPoint) const;
        Vector ReturnStressAtRegularFailureZone(const Vector& rPrincipalTrialStressVector,
                                                const CoulombYieldFunction& rCoulombYieldFunction,
                                                const double FrictionAngle,
                                                const double Cohesion) const ;

        double CalculateApex(const double FrictionAngle, const double Cohesion) const;
        Vector CalculateCornerPoint(const double FrictionAngle, const double Cohesion, const double TensionCutoff,
                                    const double Apex) const;

    private:
        double mStateVariable;
        std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
        Vector mStressVector;
        Vector mStressVectorFinalized;
        Vector mStrainVectorFinalized;

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
