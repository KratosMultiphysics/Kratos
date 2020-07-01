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
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, int, AUX_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Matrix, ROM_BASIS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, double, HROM_WEIGHT )

    // Modal derivative variables
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, unsigned int, BUILD_LEVEL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, Vector, EIGENVALUE_VECTOR)
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, std::size_t, BASIS_I )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, std::size_t, BASIS_J )
    KRATOS_DEFINE_APPLICATION_VARIABLE( ROM_APPLICATION, std::size_t, DERIVATIVE_INDEX )

}

#endif	/* KRATOS_ROM_APPLICATION_VARIABLES_H_INCLUDED */
