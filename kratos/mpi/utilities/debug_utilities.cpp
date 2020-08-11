//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Carlos A. Roig
//                   Riccardo Rossi
//


#include "mpi/utilities/debug_utilities.h"

namespace Kratos {

template<class TVariableType>
bool MpiDebugUtilities::InteralCmpEq(const TVariableType& rVar1, const TVariableType& rVar2) {
    return rVar1 == rVar2;
}

} // namespace Kratos
