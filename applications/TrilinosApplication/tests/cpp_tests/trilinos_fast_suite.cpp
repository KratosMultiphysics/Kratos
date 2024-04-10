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

#include "trilinos_fast_suite.h"

namespace Kratos::Testing {
KratosTrilinosFastSuite::KratosTrilinosFastSuite()
    : KratosMPICoreFastSuite() {
  if (!this->mKernel.IsImported("TrilinosApplication")) {
    mpTrilinosApp = std::make_shared<KratosTrilinosApplication>();
    this->mKernel.ImportApplication(mpTrilinosApp);
  }
}

} // namespace Kratos::Testing

// Create a custom main with the MPI environment and custom listeners for the test output
int main(int argc, char* argv[]) 
{
    return Kratos::Testing::GTestMain::InitializeMPIKernel(argc, argv);
}