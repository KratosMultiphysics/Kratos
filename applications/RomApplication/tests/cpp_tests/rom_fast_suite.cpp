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

#include "rom_fast_suite.h"

namespace Kratos::Testing {
KratosRomFastSuite::KratosRomFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("RomApplication")) {
    mpRomApp = std::make_shared<KratosRomApplication>();
    this->mKernel.ImportApplication(mpRomApp);
  }
}

} // namespace Kratos::Testing
