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
    Pointer       Clone() const override;
    SizeType      WorkingSpaceDimension() override;
    SizeType      GetStrainSize() const override;
    StressMeasure GetStressMeasure() override;
    bool          IsIncremental() override;
    int           Check(const Properties&   rMaterialProperties,
                        const GeometryType& rElementGeometry,
                        const ProcessInfo&  rCurrentProcessInfo) const override;
    void          CalculateMaterialResponseCauchy(Parameters& rValues) override;
    bool          RequiresInitializeMaterialResponse() override;
    void          FinalizeMaterialResponseCauchy(Parameters& rValues) override;
    void InitializeMaterial(const Properties&, const GeometryType&, const Vector&) override;

private:
    Vector mPreviousRelativeDisplacement;
    Vector mPreviousTraction;
};

} // namespace Kratos
