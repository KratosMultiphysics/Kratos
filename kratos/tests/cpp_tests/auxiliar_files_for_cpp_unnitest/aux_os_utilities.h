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
//

// Project includes
#include "utilities/os_utilities.h"

namespace Kratos {
namespace Testing {

static inline std::string AuxiliarGetCurrentWorkingDir()
{
    return OSUtilities::GetCurrentWorkingDir();
}

}   // namespace Testing
}  // namespace Kratos.
