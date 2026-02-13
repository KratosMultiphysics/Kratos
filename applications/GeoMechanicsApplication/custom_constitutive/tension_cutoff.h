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

namespace Kratos
{

namespace Geo
{
class PrincipalStresses;
class SigmaTau;
} // namespace Geo

class KRATOS_API(GEO_MECHANICS_APPLICATION) TensionCutoff : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TensionCutoff);

    TensionCutoff() = default;

    explicit TensionCutoff(double TensileStrength);

    [[nodiscard]] double GetTensileStrength() const;

    [[nodiscard]] double YieldFunctionValue(const Vector& rSigmaTau) const override;
    [[nodiscard]] double YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const;
    [[nodiscard]] double YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::SigmaTau&) const;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    double mTensileStrength = 0.0;

}; // Class TensionCutoff

} // namespace Kratos