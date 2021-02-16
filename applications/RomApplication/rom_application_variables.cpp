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
    KRATOS_CREATE_VARIABLE( double, HROM_WEIGHT )
    KRATOS_CREATE_VARIABLE( TMapPhiType, MAP_PHI )

    // Modal derivative variables
    KRATOS_CREATE_VARIABLE(unsigned int, BUILD_LEVEL)
    KRATOS_CREATE_VARIABLE(Vector, EIGENVALUE_VECTOR)
    KRATOS_CREATE_VARIABLE(std::size_t, BASIS_I)
    KRATOS_CREATE_VARIABLE(std::size_t, BASIS_J)
    KRATOS_CREATE_VARIABLE(std::size_t, DERIVATIVE_INDEX )
    KRATOS_CREATE_VARIABLE(double, MODAL_COORDINATE )
}
