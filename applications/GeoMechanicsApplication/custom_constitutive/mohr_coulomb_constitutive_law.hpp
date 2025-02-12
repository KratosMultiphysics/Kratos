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
#include "includes/constitutive_law.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "includes/serializer.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) MohrCoulombConstitutiveLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MohrCoulombConstitutiveLaw);

    MohrCoulombConstitutiveLaw();
    //explicit MohrCoulombConstitutiveLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    ~MohrCoulombConstitutiveLaw() override = default;

    ConstitutiveLaw::Pointer Clone() const override;


    using ConstitutiveLaw::GetValue;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    using ConstitutiveLaw::SetValue;
    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMohrCoulomb(const Properties& rProp, Vector& rCautchyStressVector);

    int FindRegionIndex(double fme, double fte);

    void CalculatePK2Stress(const Vector& rStrainVector, Vector& rStressVector, ConstitutiveLaw::Parameters& rValues);
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues);
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues);

    // Member Variables
    double mStateVariable;

private:
    //std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector mStressVector;
    Vector mStrainVectorFinalized;
    Vector mStressVectorFinalized;
    Vector mDeltaStrainVector;

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