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
    CoulombWithTensionCutOffImpl(double FrictionAngleInRadians,
                                 double Cohesion,
                                 double DilatancyAngleInRadians,
                                 double TensileStrength);

    [[nodiscard]] bool IsAdmissibleSigmaTau(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] Vector DoReturnMapping(const Properties& rProperties, const Vector& rTrialSigmaTau) const;

private:
    CoulombYieldSurface mCoulombYieldSurface;
    TensionCutoff       mTensionCutOff;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);
};

} // namespace Kratos
