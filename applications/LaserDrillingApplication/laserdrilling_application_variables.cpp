//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#include "laserdrilling_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, NO2)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(TRANSPORT_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_VELOCITAT_VENT)
KRATOS_CREATE_VARIABLE(bool, IS_POTENTIAL_FLOW_STEP)
KRATOS_CREATE_VARIABLE(double, VELOCITAT_VENT)
KRATOS_CREATE_VARIABLE(double, DIRECCIO_VENT)
KRATOS_CREATE_VARIABLE(double, METEO_DIRECCIO_VENT)
KRATOS_CREATE_VARIABLE(std::string, STARTING_DATE)
KRATOS_CREATE_VARIABLE(int, SIMULATION_DURATION_IN_DAYS)
KRATOS_CREATE_VARIABLE(bool, WIND_AUTOMATIC_PROCESS)
KRATOS_CREATE_VARIABLE(bool, POLLUTANT_AUTOMATIC_PROCESS)
KRATOS_CREATE_VARIABLE(std::string, CITY)
KRATOS_CREATE_VARIABLE(int, CASE_ID)
KRATOS_CREATE_VARIABLE(bool, IN_PRODUCTION)

}
