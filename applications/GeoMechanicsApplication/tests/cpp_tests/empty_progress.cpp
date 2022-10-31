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

#include "empty_progress.h"

namespace empty_progress
{
    void emptyProgress(double progress) {}
    void emptyLog(char* log) {}
    bool emptyCancel() {
        return false;
    }
}