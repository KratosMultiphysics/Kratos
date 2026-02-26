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

#include "custom_constitutive/principal_stresses.hpp"
#include "includes/kratos_export_api.h"
#include "includes/smart_pointers.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

namespace Geo
{
class PrincipalStresses;
class SigmaTau;
} // namespace Geo

class KRATOS_API(GEO_MECHANICS_APPLICATION) TensionCutoff
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TensionCutoff);

    TensionCutoff() = default;

    explicit TensionCutoff(double TensileStrength);

    [[nodiscard]] double GetTensileStrength() const;

    [[nodiscard]] double YieldFunctionValue(const Geo::SigmaTau& rSigmaTau) const;
    [[nodiscard]] double YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::SigmaTau&,
                                                  Geo::PrincipalStresses::AveragingType AveragingType =
                                                      Geo::PrincipalStresses::AveragingType::NO_AVERAGING) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::PrincipalStresses&,
                                                  Geo::PrincipalStresses::AveragingType AveragingType =
                                                      Geo::PrincipalStresses::AveragingType::NO_AVERAGING) const;

    [[nodiscard]] double CalculatePlasticMultiplier(const Geo::SigmaTau& rTrialSigmaTau,
                                                    const Vector&        rDerivativeOfFlowFunction,
                                                    const Matrix&        rElasticMatrix) const;
    [[nodiscard]] double CalculatePlasticMultiplier(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                    const Vector& rDerivativeOfFlowFunction,
                                                    const Matrix& rElasticMatrix) const;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);

    double mTensileStrength = 0.0;

}; // Class TensionCutoff

} // namespace Kratos