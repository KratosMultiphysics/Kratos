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

#include "mappings_fast_suite.h"

namespace Kratos::Testing {
KratosMappingFastSuite::KratosMappingFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("MappingApplication")) {
    mpMappingApp = std::make_shared<KratosMappingApplication>();
    this->mKernel.ImportApplication(mpMappingApp);
  }
}

} // namespace Kratos::Testing
