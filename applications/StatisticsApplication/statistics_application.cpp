//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

// System includes

// External includes

// Project includes
#include "statistics_application.h"
#include "statistics_application_variables.h"

namespace Kratos
{
KratosStatisticsApplication::KratosStatisticsApplication()
    : KratosApplication("StatisticsApplication")
{
}

void KratosStatisticsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosStatisticsApplication..." << std::endl;

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_SUM)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_MEAN)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_VARIANCE)

    KRATOS_REGISTER_VARIABLE(VELOCITY_NORM)
    KRATOS_REGISTER_VARIABLE(PRESSURE_NORM)
    KRATOS_REGISTER_VARIABLE(PRESSURE_SUM)
    KRATOS_REGISTER_VARIABLE(PRESSURE_MEAN)
    KRATOS_REGISTER_VARIABLE(PRESSURE_VARIANCE)
}
} // namespace Kratos.
