//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
//

// Project includes
#include "trilinos_fast_suite.h"

namespace Kratos::Testing {

KratosTrilinosApplicationMPITestSuite::KratosTrilinosApplicationMPITestSuite()
    : KratosMPICoreFastSuite() {
  mpTrilinosApp = std::make_shared<KratosTrilinosApplication>();
  this->ImportApplicationIntoKernel(mpTrilinosApp);
}

} // namespace Kratos::Testing