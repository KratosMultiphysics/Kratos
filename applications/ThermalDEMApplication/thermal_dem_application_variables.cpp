//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
//

// Project includes
#include "thermal_dem_application_variables.h"

namespace Kratos
{
  // Variables
  KRATOS_CREATE_VARIABLE(ThermalDEMIntegrationScheme::Pointer, THERMAL_INTEGRATION_SCHEME_POINTER)
  KRATOS_CREATE_VARIABLE(NumericalIntegrationMethod::Pointer,  NUMERICAL_INTEGRATION_METHOD_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       DIRECT_CONDUCTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       INDIRECT_CONDUCTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       CONVECTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatExchangeMechanism::Pointer,       RADIATION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(HeatGenerationMechanism::Pointer,     FRICTION_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(RealContactModel::Pointer,            REAL_CONTACT_MODEL_POINTER)
  KRATOS_CREATE_VARIABLE(std::string,                          THERMAL_INTEGRATION_SCHEME_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          NUMERICAL_INTEGRATION_METHOD_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          DIRECT_CONDUCTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          INDIRECT_CONDUCTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          CONVECTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          RADIATION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          FRICTION_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          REAL_CONTACT_MODEL_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          VORONOI_METHOD_NAME)
  KRATOS_CREATE_VARIABLE(std::string,                          POROSITY_METHOD_NAME)
  KRATOS_CREATE_VARIABLE(bool,                                 MOTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 DIRECT_CONDUCTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 INDIRECT_CONDUCTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 CONVECTION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 RADIATION_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 FRICTION_HEAT_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 REAL_CONTACT_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 TEMPERATURE_DEPENDENT_RADIUS_OPTION)
  KRATOS_CREATE_VARIABLE(bool,                                 FIXED_TEMPERATURE)
  KRATOS_CREATE_VARIABLE(bool,                                 ADIABATIC)
  KRATOS_CREATE_VARIABLE(int,                                  THERMAL_FREQUENCY)
  //KRATOS_CREATE_VARIABLE(double,                             HEATFLUX)             // defined in DEMApp
  //KRATOS_CREATE_VARIABLE(double,                             THERMAL_CONDUCTIVITY) // defined in DEMApp
  KRATOS_CREATE_VARIABLE(double,                               REAL_YOUNG_MODULUS_RATIO)
  KRATOS_CREATE_VARIABLE(double,                               HEATSOURCE)
  KRATOS_CREATE_VARIABLE(double,                               MIN_CONDUCTION_DISTANCE)
  KRATOS_CREATE_VARIABLE(double,                               MAX_CONDUCTION_DISTANCE)
  KRATOS_CREATE_VARIABLE(double,                               ISOTHERMAL_CORE_RADIUS)
  KRATOS_CREATE_VARIABLE(double,                               MAX_RADIATION_DISTANCE)
  KRATOS_CREATE_VARIABLE(double,                               FRICTION_HEAT_CONVERSION)
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

  // Flags (starting at 1000 to avoid conflicts with DEMFlags)
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, IS_ADIABATIC,                     1000);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_MOTION,                       1001);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_FIXED_TEMPERATURE,            1002);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_DIRECT_CONDUCTION,            1003);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_INDIRECT_CONDUCTION,          1004);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_CONVECTION,                   1005);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_RADIATION,                    1006);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_FRICTION_HEAT,                1007);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_REAL_CONTACT,                 1008);
  KRATOS_CREATE_LOCAL_FLAG(DEMThermalFlags, HAS_TEMPERATURE_DEPENDENT_RADIUS, 1009);

} // namespace Kratos
