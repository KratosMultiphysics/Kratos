// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#pragma once

namespace Kratos::MicroClimateConstants
{

constexpr static auto AirDensity = 1.18;
constexpr static auto AirHeatCapacity = 1004.67;
constexpr static auto LatentEvaporationHeat = 2.45e6;
constexpr static auto PsychometricConstant = 0.63;
constexpr static auto EffectiveEmissivity = 0.95;
constexpr static auto SurfaceResistance = 30.0;
constexpr static auto MeasurementHeight = 10.0;

constexpr static auto RoughnessLayerResistance = 30.0;
constexpr static auto RoughnessLayerHeight = 10.0;
constexpr static auto RoughnessHeight = 1.0;

constexpr static auto VonNeumannCoefficient = 0.4;
constexpr static auto BoltzmannCoefficient = 5.67e-8;
constexpr static auto GravitationalAcceleration = 9.81;

} // namespace Kratos::MicroClimateConstants