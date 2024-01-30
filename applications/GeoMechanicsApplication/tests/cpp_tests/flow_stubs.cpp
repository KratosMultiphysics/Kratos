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
    void emptyProgress(double progress) {}
    void emptyLog(const char* log) {}
    bool emptyCancel() {
        return false;
    }
}