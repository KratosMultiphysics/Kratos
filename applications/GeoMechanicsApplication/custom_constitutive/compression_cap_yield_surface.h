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
//                   Anne van de Graaf
//

#pragma once

#include "custom_constitutive/yield_surface.h"
#include "geo_aliases.h"
#include "includes/properties.h"

namespace Kratos
{

class CheckProperties;

class KRATOS_API(GEO_MECHANICS_APPLICATION) CompressionCapYieldSurface : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CompressionCapYieldSurface);

    CompressionCapYieldSurface();
    explicit CompressionCapYieldSurface(const Properties& rMaterialProperties);

    [[nodiscard]] double GetCapSize() const;
    [[nodiscard]] double GetPreconsolidationStress() const;

    [[nodiscard]] double YieldFunctionValue(const Vector& rSigmaTau) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&) const override;

private:
    void InitializeKappaDependentFunctions();
    void CheckMaterialProperties() const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    double                      mKappa = 0.0;
    Properties                  mMaterialProperties;
    Geo::KappaDependentFunction mCapSizeCalculator;
    Geo::KappaDependentFunction mPreconsolidationStressCalculator;

}; // Class CompressionCapYieldSurface

} // namespace Kratos