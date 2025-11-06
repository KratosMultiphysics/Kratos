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

class KRATOS_API(GEO_MECHANICS_APPLICATION) CoulombYieldSurface : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CoulombYieldSurface);

    using KappaDependentFunction = std::function<double(double)>;

    enum class CoulombAveragingType {
        NO_AVERAGING,
        LOWEST_PRINCIPAL_STRESSES,
        HIGHEST_PRINCIPAL_STRESSES
    };

    CoulombYieldSurface();
    explicit CoulombYieldSurface(Properties MaterialProperties);

    [[nodiscard]] double GetFrictionAngleInRadians() const;
    [[nodiscard]] double GetCohesion() const;
    [[nodiscard]] double GetDilatancyAngleInRadians() const;
    [[nodiscard]] double GetKappa() const;

    void SetKappa(const double kappa);

    [[nodiscard]] double YieldFunctionValue(const Vector& rSigmaTau) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&, CoulombAveragingType AveragingType) const;

private:
    void InitializeKappaDependentFunctions();

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    double                 mKappa = 0.0;
    Properties             mMaterialProperties;
    KappaDependentFunction mFrictionAngleCalculator;
    KappaDependentFunction mCohesionCalculator;
    KappaDependentFunction mDilatancyAngleCalculator;

}; // Class CoulombYieldSurface

} // namespace Kratos