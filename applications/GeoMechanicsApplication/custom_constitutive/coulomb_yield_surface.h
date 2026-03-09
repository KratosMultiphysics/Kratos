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

#include "custom_constitutive/principal_stresses.hpp"
#include "custom_constitutive/sigma_tau.hpp"
#include "geo_aliases.h"
#include "includes/properties.h"

#include <functional>

namespace Kratos
{

class CheckProperties;

namespace Geo
{
class SigmaTau;
class PQ;
} // namespace Geo

class KRATOS_API(GEO_MECHANICS_APPLICATION) CoulombYieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CoulombYieldSurface);

    CoulombYieldSurface();
    explicit CoulombYieldSurface(const Properties& rMaterialProperties);

    [[nodiscard]] double GetFrictionAngleInRadians() const;
    [[nodiscard]] double GetCohesion() const;
    [[nodiscard]] double GetDilatancyAngleInRadians() const;
    [[nodiscard]] double GetKappa() const;
    void                 SetKappa(double kappa);

    [[nodiscard]] double YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const;
    [[nodiscard]] double YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const;
    [[nodiscard]] double YieldFunctionValue(const Geo::PQ& rPQ) const;

    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::SigmaTau&,
                                                  Geo::PrincipalStresses::AveragingType AveragingType =
                                                      Geo::PrincipalStresses::AveragingType::NO_AVERAGING) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::PrincipalStresses&,
                                                  Geo::PrincipalStresses::AveragingType AveragingType =
                                                      Geo::PrincipalStresses::AveragingType::NO_AVERAGING) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::PQ& rPQ,
                                                  Geo::PrincipalStresses::AveragingType AveragingType) const;

    [[nodiscard]] Geo::SigmaTau CalculateApex() const;

    [[nodiscard]] double CalculatePlasticMultiplier(const Geo::SigmaTau& rTrialSigmaTau,
                                                    const Vector&        rDerivativeOfFlowFunction,
                                                    const Matrix&        rElasticMatrix) const;
    [[nodiscard]] double CalculatePlasticMultiplier(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                    const Vector& rDerivativeOfFlowFunction,
                                                    const Matrix& rElasticMatrix) const;
    [[nodiscard]] double CalculateEquivalentPlasticStrainIncrement(const Geo::SigmaTau& rTrialSigmaTau,
                                                                   const Matrix& rElasticMatrix,
                                                                   Geo::PrincipalStresses::AveragingType AveragingType) const;
    [[nodiscard]] double CalculateEquivalentPlasticStrainIncrement(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                                   const Matrix& rElasticMatrix,
                                                                   Geo::PrincipalStresses::AveragingType AveragingType) const;

private:
    void InitializeKappaDependentFunctions();
    void CheckMaterialProperties() const;
    void CheckHardeningCoefficients(const Variable<Vector>& rCoefficientsVariable,
                                    const CheckProperties&  rChecker) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);

    double                      mKappa = 0.0;
    Properties                  mMaterialProperties;
    Geo::KappaDependentFunction mFrictionAngleCalculator;
    Geo::KappaDependentFunction mCohesionCalculator;
    Geo::KappaDependentFunction mDilatancyAngleCalculator;
}; // Class CoulombYieldSurface

} // namespace Kratos