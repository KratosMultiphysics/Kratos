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

#include "fsi_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosFSIFastSuite : public KratosCoreFastSuite {
public:
  KratosFSIFastSuite();

private:
  KratosFSIApplication::Pointer mpFSIApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
