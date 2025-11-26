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
    explicit CoulombYieldSurface(const Properties& rMaterialProperties);

    [[nodiscard]] double GetFrictionAngleInRadians() const;
    [[nodiscard]] double GetCohesion() const;
    [[nodiscard]] double GetDilatancyAngleInRadians() const;
    [[nodiscard]] double GetKappa() const;
    void                 SetKappa(double kappa);

    [[nodiscard]] double YieldFunctionValue(const Vector& rSigmaTau) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&, CoulombAveragingType AveragingType) const;

    [[nodiscard]] double CalculateApex() const;
    [[nodiscard]] double CalculatePlasticMultiplier(const Vector& rSigmaTau,
                                                    const Vector& rDerivativeOfFlowFunction) const;
    [[nodiscard]] double CalculateEquivalentPlasticStrainIncrement(const Vector& rSigmaTau,
                                                                   CoulombAveragingType AveragingType) const;

private:
    void InitializeKappaDependentFunctions();
    void CheckMaterialProperties();

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