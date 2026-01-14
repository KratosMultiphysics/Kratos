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

#include "custom_constitutive/yield_surface.h"
#include "includes/properties.h"

#include <functional>

namespace Kratos
{

class CheckProperties;

class KRATOS_API(GEO_MECHANICS_APPLICATION) CompressionCapYieldSurface : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CompressionCapYieldSurface);

    using KappaDependentFunction = std::function<double(double)>;

    CompressionCapYieldSurface();
    explicit CompressionCapYieldSurface(const Properties& rMaterialProperties);

    [[nodiscard]] double GetCapSize() const;
    [[nodiscard]] double GetCapLocation() const;
    [[nodiscard]] double GetKappa() const;
    void                 SetKappa(double kappa);

    [[nodiscard]] double YieldFunctionValue(const Vector& rSigmaTau) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&) const override;

private:
    void   InitializeKappaDependentFunctions();
    void   CheckMaterialProperties() const;
    double GetCapSize(const Properties& rProperties);
    double GetCapLocation(const Properties& rProperties);

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    double                 mKappa = 0.0;
    Properties             mMaterialProperties;
    KappaDependentFunction mCapSizeCalculator;
    KappaDependentFunction mCapLocationCalculator;

}; // Class CompressionCapYieldSurface

} // namespace Kratos