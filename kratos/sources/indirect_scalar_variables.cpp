//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/variables.h"

// Include base h
#include "includes/indirect_scalar_variables.h"

namespace Kratos
{
    // Zero scalar variable
    Kratos::IndirectScalarVariable INDIRECT_SCALAR_ZERO;

    // Adjoint variables variable
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_1_X)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_1_Y)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_1_Z)


    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_2_X)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_2_Y)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_2_Z)

    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_3_X)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_3_Y)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_VECTOR_3_Z)

    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(AUX_ADJOINT_VECTOR_1_X)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(AUX_ADJOINT_VECTOR_1_Y)
    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(AUX_ADJOINT_VECTOR_1_Z)

    KRATOS_CREATE_INDIRECT_SCALAR_VARIABLE(ADJOINT_SCALAR_1)

}  // namespace Kratos.

