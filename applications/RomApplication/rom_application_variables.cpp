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
//                  Altug Emiroglu, https://github.com/emiroglu
//
//


#include "rom_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE( int, AUX_ID )
    KRATOS_CREATE_VARIABLE( Matrix, ROM_BASIS )

    // Modal derivative variables
    KRATOS_CREATE_VARIABLE(int, BUILD_LEVEL)
    KRATOS_CREATE_VARIABLE(Vector, EIGENVALUE_VECTOR)
    KRATOS_CREATE_VARIABLE(int, EIGENVALUE_I)
    KRATOS_CREATE_VARIABLE(int, EIGENVALUE_J)
}
