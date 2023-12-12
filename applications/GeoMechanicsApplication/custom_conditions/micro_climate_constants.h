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

constexpr static double air_density = 1.18;
constexpr static double air_heat_capacity = 1004.67;
constexpr static double roughness_layer_resistance = 30.0;
constexpr static double latent_evaporation_heat = 2.45e6;
constexpr static double water_density = 1e3;
constexpr static double psychometric_constant = 0.63;
constexpr static double surface_resistance = 30.0;

constexpr static double roughness_layer_height = 10.0;
constexpr static double von_neuman_coefficient = 0.4;
constexpr static double measurement_height = 10.0;
constexpr static double roughness_height = 1.0;
constexpr static double gravitational_acceleration = 9.81;

constexpr auto effective_emissivity = 0.95;
constexpr auto boltzmann_coefficient = 5.67e-8;

} // namespace Kratos::MicroClimateConstants