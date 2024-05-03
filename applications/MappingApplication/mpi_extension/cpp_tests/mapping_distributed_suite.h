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

#pragma once

#include "mapping_application.h"
#include "mpi/testing/mpi_testing.h"

namespace Kratos::Testing {

class KratosMappingMPIFastSuite : public KratosMPICoreFastSuite {
public:
  KratosMappingMPIFastSuite();

private:
  KratosMappingApplication::Pointer mpMappingApp;
};

} // namespace Kratos::Testing
