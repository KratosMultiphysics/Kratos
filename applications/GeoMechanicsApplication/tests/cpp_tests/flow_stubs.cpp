// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Carlos Lubbers
//

#include "flow_stubs.h"

namespace flow_stubs
{
void emptyProgress([[maybe_unused]] double progress)
{
    // deliberately empty as the name says
}

void emptyLog([[maybe_unused]] const char* log)
{
    // deliberately empty as the name says
}

bool emptyCancel() { return false; }
} // namespace flow_stubs