//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//
//


#include "rom_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE( Matrix, ROM_BASIS )
    KRATOS_CREATE_VARIABLE( Matrix, ROM_LEFT_BASIS )
    KRATOS_CREATE_VARIABLE( double, HROM_WEIGHT )
    KRATOS_CREATE_VARIABLE( Vector, ROM_SOLUTION_INCREMENT )
}
