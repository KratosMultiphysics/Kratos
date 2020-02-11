//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/kratos_filesystem.h"
#include "utilities/os_utilities.h"

namespace Kratos
{
namespace OSUtilities
{
std::string GetCurrentWorkingDir()
{
    KRATOS_TRY

    return filesystem::current_path().string();

    KRATOS_CATCH("")
}

} // namespace OSUtilities
} // namespace Kratos
