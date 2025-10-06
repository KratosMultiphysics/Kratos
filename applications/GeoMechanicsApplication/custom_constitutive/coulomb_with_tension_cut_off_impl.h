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

    [[nodiscard]] bool   IsAdmissibleSigmaTau(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] Vector DoReturnMapping(const Properties& rProperties,
                                         const Vector&     rTrialSigmaTau,
                                         CoulombYieldSurface::CoulombAveragingType AveragingType);

private:
    CoulombYieldSurface mCoulombYieldSurface;
    TensionCutoff       mTensionCutOff;
    double mEquivalentPlasticStrain = 0.0;

    double UpdateFrictionAngle(const Properties& rProperties, double kappa) const;
    double UpdateCohesion(const Properties& rProperties, double kappa) const;
    double UpdateDilatancyAngle(const Properties& rProperties, double kappa) const;
    double CalculateEquivalentPlasticStrain(const Vector& rSigmaTau,
        CoulombYieldSurface::CoulombAveragingType AveragingType,
        double lambda) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);
};

} // namespace Kratos
