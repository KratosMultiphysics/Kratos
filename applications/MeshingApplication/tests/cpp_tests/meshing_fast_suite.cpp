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

#include "meshing_fast_suite.h"

namespace Kratos::Testing {
KratosMeshingFastSuite::KratosMeshingFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("MeshingApplication")) {
    mpMeshingApp = std::make_shared<KratosMeshingApplication>();
    this->mKernel.ImportApplication(mpMeshingApp);
  }
}

} // namespace Kratos::Testing
