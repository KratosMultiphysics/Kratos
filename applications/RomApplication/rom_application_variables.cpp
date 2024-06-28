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
    KRATOS_CREATE_VARIABLE( Vector, ROM_SOLUTION_BASE )
    KRATOS_CREATE_VARIABLE( Vector, ROM_SOLUTION_TOTAL )
    KRATOS_CREATE_VARIABLE( Vector, SOLUTION_BASE )
    KRATOS_CREATE_VARIABLE( vector<Matrix>, SVD_PHI_MATRICES )
    KRATOS_CREATE_VARIABLE( vector<Matrix>, NN_LAYERS )
    KRATOS_CREATE_VARIABLE( Vector, SOLUTION_REFERENCE )
    KRATOS_CREATE_VARIABLE( bool, UPDATE_PHI_EFFECTIVE_BOOL )
    
    // KRATOS_CREATE_VARIABLE( vector<EigenDynamicMatrix>, SVD_PHI_MATRICES )
    // KRATOS_CREATE_VARIABLE( vector<EigenDynamicMatrix>, NN_LAYERS )
    // KRATOS_CREATE_VARIABLE( EigenDynamicVector, SOLUTION_REFERENCE )
    // KRATOS_CREATE_VARIABLE( bool, UPDATE_PHI_EFFECTIVE_BOOL )
}
