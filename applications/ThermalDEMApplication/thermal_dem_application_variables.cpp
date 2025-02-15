//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// Project includes
#include "thermal_dem_application_variables.h"

namespace Kratos
{
  // Variables
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       DIRECT_CONDUCTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       INDIRECT_CONDUCTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       CONVECTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       RADIATION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatGenerationMechanism::Pointer,     GENERATION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(RealContactModel::Pointer,            REAL_CONTACT_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(ThermalDEMIntegrationScheme::Pointer, THERMAL_INTEGRATION_SCHEME_POINTER)
  KRATOS_CREATE_VARIABLE(NumericalIntegrationMethod::Pointer,  NUMERICAL_INTEGRATION_METHOD_POINTER)
  KRATOS_CREATE_VARIABLE(std::string,                          DIRECT_CONDUCTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          INDIRECT_CONDUCTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          CONVECTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          RADIATION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          REAL_CONTACT_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          THERMAL_INTEGRATION_SCHEME_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          NUMERICAL_INTEGRATION_METHOD_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          VORONOI_METHOD_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          POROSITY_METHOD_NAME)
  KRATOS_CREATE_VARIABLE(bool,                                 AUTO_SOLVE_FREQUENCY_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 COMPUTE_FORCES_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 COMPUTE_MOTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 DIRECT_CONDUCTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 INDIRECT_CONDUCTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 CONVECTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 RADIATION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 HEAT_GENERATION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 GENERATION_SLIDING_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 GENERATION_ROLLING_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 GENERATION_DAMPING_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 HEAT_MAP_GENERATION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 REAL_CONTACT_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 FIXED_TEMPERATURE)
  KRATOS_CREATE_VARIABLE(bool,                                 ADIABATIC)
  KRATOS_CREATE_VARIABLE(int,                                  THERMAL_FREQUENCY)
  KRATOS_CREATE_VARIABLE(double,                               HEATFLUX)
  KRATOS_CREATE_VARIABLE(double,                               THERMAL_CONDUCTIVITY)
  KRATOS_CREATE_VARIABLE(double,                               DEFORMATION_RATE)
  KRATOS_CREATE_VARIABLE(double,                               DEFORMATION_RATE_START_TIME)
  KRATOS_CREATE_VARIABLE(double,                               DEFORMATION_RATE_STOP_TIME)
  KRATOS_CREATE_VARIABLE(double,                               REAL_YOUNG_MODULUS_RATIO)
  KRATOS_CREATE_VARIABLE(double,                               HEATSOURCE)
  KRATOS_CREATE_VARIABLE(double,                               MIN_CONDUCTION_DISTANCE)
  KRATOS_CREATE_VARIABLE(double,                               MAX_CONDUCTION_DISTANCE)
  KRATOS_CREATE_VARIABLE(double,                               CONDUCTION_RADIUS)
  KRATOS_CREATE_VARIABLE(double,                               ISOTHERMAL_CORE_RADIUS)
  KRATOS_CREATE_VARIABLE(double,                               MAX_RADIATION_DISTANCE)
  KRATOS_CREATE_VARIABLE(double,                               HEAT_GENERATION_RATIO)
  KRATOS_CREATE_VARIABLE(double,                               AVERAGE_POROSITY)
  KRATOS_CREATE_VARIABLE(double,                               ALPHA_SHAPE_PARAMETER)
  KRATOS_CREATE_VARIABLE(double,                               INTEGRAL_TOLERANCE)
  KRATOS_CREATE_VARIABLE(double,                               FLUID_LAYER_THICKNESS)
  KRATOS_CREATE_VARIABLE(double,                               FLUID_DENSITY)
  KRATOS_CREATE_VARIABLE(double,                               FLUID_VISCOSITY)
  KRATOS_CREATE_VARIABLE(double,                               FLUID_THERMAL_CONDUCTIVITY)
  KRATOS_CREATE_VARIABLE(double,                               FLUID_HEAT_CAPACITY)
  KRATOS_CREATE_VARIABLE(double,                               FLUID_TEMPERATURE)
  KRATOS_CREATE_VARIABLE(Vector,                               FLUID_VELOCITY)
  KRATOS_CREATE_VARIABLE(Vector,                               HEAT_MAP_COORDINATES_1)
  KRATOS_CREATE_VARIABLE(Vector,                               HEAT_MAP_COORDINATES_2)
  KRATOS_CREATE_VARIABLE(Vector,                               HEAT_MAP_SUBDIVISIONS)

  // Flags (starting from the last value of the DEM Appication)
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, IS_ADIABATIC, 18);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, IS_SINTERING, 19);

} // namespace Kratos
