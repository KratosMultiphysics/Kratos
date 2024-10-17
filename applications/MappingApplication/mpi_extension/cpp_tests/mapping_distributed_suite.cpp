//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//

// External includes
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// Project includes
#include "mapping_distributed_suite.h"

namespace Kratos::Testing {

KratosMappingMPIFastSuite::KratosMappingMPIFastSuite()
    : KratosMPICoreFastSuite() {
  mpMappingApp = std::make_shared<KratosMappingApplication>();
  this->ImportApplicationIntoKernel(mpMappingApp);
}

} // namespace Kratos::Testing