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
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "coulomb_yield_surface.h"
#include "tension_cutoff.h"

namespace Kratos
{

class Properties;
class Serializer;

class CoulombWithTensionCutOffImpl
{
public:
    CoulombWithTensionCutOffImpl() = default;
    explicit CoulombWithTensionCutOffImpl(const Properties& rMaterialProperties);

    [[nodiscard]] bool   IsAdmissibleSigmaTau(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] Vector DoReturnMapping(const Vector& rTrialSigmaTau,
                                         CoulombYieldSurface::CoulombAveragingType AveragingType);
private:
    CoulombYieldSurface mCoulombYieldSurface;
    TensionCutoff       mTensionCutOff;

    double CalculateEquivalentPlasticStrain(const Vector&                             rSigmaTau,
                                            CoulombYieldSurface::CoulombAveragingType AveragingType,
                                            double                                    lambda) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);
};

} // namespace Kratos
