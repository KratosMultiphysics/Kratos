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

KRATOS_CREATE_VARIABLE(double, THERMAL_ENERGY)
KRATOS_CREATE_VARIABLE(double, THERMAL_ENERGY_PER_VOLUME)
KRATOS_CREATE_VARIABLE(double, THERMAL_DECOMPOSITION)
KRATOS_CREATE_VARIABLE(double, DECOMPOSITION_LAW_REFERENCE_TEMPERATURE)
KRATOS_CREATE_VARIABLE(double, DECOMPOSITION_LAW_CONSTANT_1)
KRATOS_CREATE_VARIABLE(double, DECOMPOSITION_LAW_CONSTANT_2)
KRATOS_CREATE_VARIABLE(double, THERMAL_ENERGY_PER_VOLUME_THRESHOLD)
KRATOS_CREATE_VARIABLE(double, ALPHA_THRESHOLD)
KRATOS_CREATE_VARIABLE(int, THERMAL_COUNTER)
KRATOS_CREATE_VARIABLE(std::string, DECOMPOSITION_LAW)

}
