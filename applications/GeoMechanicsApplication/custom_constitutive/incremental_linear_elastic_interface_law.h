// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "includes/constitutive_law.h"
#include "includes/kratos_export_api.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoIncrementalLinearElasticInterfaceLaw : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;

    [[nodiscard]] Pointer  Clone() const override;
    SizeType               WorkingSpaceDimension() override;
    [[nodiscard]] SizeType GetStrainSize() const override;
    Vector&                GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    using BaseType::GetValue;
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;
    using BaseType::CalculateValue;
    StressMeasure GetStressMeasure() override;
    bool          IsIncremental() override;
    void InitializeMaterial(const Properties&, const GeometryType&, const Vector&) override;
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;
    bool RequiresInitializeMaterialResponse() override;
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;
    int  Check(const Properties&   rMaterialProperties,
               const GeometryType& rElementGeometry,
               const ProcessInfo&  rCurrentProcessInfo) const override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    Vector mPreviousRelativeDisplacement;
    Vector mPreviousTraction;
};

} // namespace Kratos
