//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//                   Simon Wenczowski
//

#include "includes/variables.h"
#include "hello.h"

namespace Kratos {

void Hello::Greet()
{
    KRATOS_INFO("Hello") << "A first greeting from the GeodataProcessingApplication: Hello World!" << std::endl;
}

}