// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include "includes/constitutive_law.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoLinearElasticLaw : public ConstitutiveLaw
{
public:
    bool RequiresInitializeMaterialResponse() override;
    bool RequiresFinalizeMaterialResponse() override;

    StrainMeasure GetStrainMeasure() override;
    StressMeasure GetStressMeasure() override;

    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>&      rThisVariable,
                           double&                      rValue) override;
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Vector>&      rThisVariable,
                           Vector&                      rValue) override;
    Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Matrix>&      rThisVariable,
                           Matrix&                      rValue) override;
    using ConstitutiveLaw::CalculateValue;

    void SetValue(const Variable<double>&, const double&, const ProcessInfo&) override;
    void SetValue(const Variable<Vector>&, const Vector&, const ProcessInfo&) override;
    using ConstitutiveLaw::SetValue;

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    void               SetConsiderDiagonalEntriesOnlyAndNoShear(bool Whether);
    [[nodiscard]] bool GetConsiderDiagonalEntriesOnlyAndNoShear() const;

protected:
    virtual void CalculateElasticMatrix(Matrix& rConstitutiveMatrix, ConstitutiveLaw::Parameters& rValues) = 0;
    virtual void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector) = 0;
    virtual void CalculatePK2Stress(const Vector&                rStrainVector,
                                    Vector&                      rStressVector,
                                    ConstitutiveLaw::Parameters& rValues) = 0;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    bool mConsiderDiagonalEntriesOnlyAndNoShear = false;
};

} // namespace Kratos
