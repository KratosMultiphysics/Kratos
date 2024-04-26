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



#if !defined(KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Matrix, ROM_BASIS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Matrix, ROM_LEFT_BASIS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, double, HROM_WEIGHT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Vector, ROM_SOLUTION_INCREMENT )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Vector, ROM_SOLUTION_BASE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Vector, ROM_SOLUTION_TOTAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Vector, SOLUTION_BASE )
}

#endif	/* KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED */
