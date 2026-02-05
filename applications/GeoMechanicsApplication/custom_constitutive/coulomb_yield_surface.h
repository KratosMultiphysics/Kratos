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

#include <functional>

namespace Kratos
{

class CheckProperties;

namespace Geo
{
class PrincipalStresses;
class SigmaTau;
} // namespace Geo

class KRATOS_API(GEO_MECHANICS_APPLICATION) CoulombYieldSurface : public YieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CoulombYieldSurface);

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
    [[nodiscard]] double YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const;
    [[nodiscard]] double YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const;

    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&) const override;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Vector&, CoulombAveragingType AveragingType) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::SigmaTau&,
                                                  CoulombAveragingType AveragingType = CoulombAveragingType::NO_AVERAGING) const;

    [[nodiscard]] double CalculateApex() const;
    [[nodiscard]] double CalculatePlasticMultiplier(const Geo::SigmaTau& rSigmaTau,
                                                    const Vector& rDerivativeOfFlowFunction) const;
    [[nodiscard]] double CalculateEquivalentPlasticStrainIncrement(const Geo::SigmaTau& rSigmaTau,
                                                                   CoulombAveragingType AveragingType) const;
    [[nodiscard]] double CalculateEquivalentPlasticStrainIncrement(const Geo::PrincipalStresses& rPrincipalStresses,
                                                                   CoulombAveragingType AveragingType) const;

private:
    void InitializeKappaDependentFunctions();
    void CheckMaterialProperties() const;
    void CheckHardeningCoefficients(const Variable<Vector>& rCoefficientsVariable,
                                    const CheckProperties&  rChecker) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    double                      mKappa = 0.0;
    Properties                  mMaterialProperties;
    Geo::KappaDependentFunction mFrictionAngleCalculator;
    Geo::KappaDependentFunction mCohesionCalculator;
    Geo::KappaDependentFunction mDilatancyAngleCalculator;
}; // Class CoulombYieldSurface

} // namespace Kratos