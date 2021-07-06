//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "integration_flags.h"

namespace Kratos
{
    KRATOS_CREATE_LOCAL_FLAG(IntegrationFlags, DO_NOT_CREATE_TESSELLATION_ON_SLAVE, 0);
    KRATOS_CREATE_LOCAL_FLAG(IntegrationFlags, USE_REVERTED_TRIM, 1);
}  // namespace Kratos.
