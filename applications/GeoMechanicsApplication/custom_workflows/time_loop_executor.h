// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include <vector>
#include <memory>

namespace Kratos
{

class Process;

class TimeLoopExecutor {
public :
    virtual ~TimeLoopExecutor() = default;
    virtual void SetProcessReferences(const std::vector<std::weak_ptr<Process>>& rProcessRefs) = 0;
};

}