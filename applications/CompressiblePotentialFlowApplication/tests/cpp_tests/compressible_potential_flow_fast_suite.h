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

#pragma once

#include "testing/testing.h"
#include "compressible_potential_flow_application.h"

namespace Kratos::Testing {

class KratosCompressiblePotentialFlowFastSuite : public KratosCoreFastSuite {
public:
  KratosCompressiblePotentialFlowFastSuite();

private:
  KratosCompressiblePotentialFlowApplication::Pointer mpCompressiblePotentialFlowApp;
};

} // namespace Kratos::Testing
