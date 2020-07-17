//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
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
    KRATOS_INFO("")
        << "    Kratos  ____  _        _   _     _   _\n"
        << "           / ___|| |_ __ _| |_(_)___| |_(_) ___ ___\n"
        << "           \\___ \\| __/ _` | __| / __| __| |/ __/ __|\n"
        << "            ___) | || (_| | |_| \\__ \\ |_| | (__\\__ \\\n"
        << "           |____/ \\__\\__,_|\\__|_|___/\\__|_|\\___|___/ "
           "Application\n"
        << "Initializing KratosStatisticsApplication..." << std::endl;

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_3D_SUM)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_3D_MEAN)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_3D_VARIANCE)

    KRATOS_REGISTER_VARIABLE(VECTOR_3D_NORM)
    KRATOS_REGISTER_VARIABLE(SCALAR_NORM)
    KRATOS_REGISTER_VARIABLE(SCALAR_SUM)
    KRATOS_REGISTER_VARIABLE(SCALAR_MEAN)
    KRATOS_REGISTER_VARIABLE(SCALAR_VARIANCE)
}
} // namespace Kratos.
