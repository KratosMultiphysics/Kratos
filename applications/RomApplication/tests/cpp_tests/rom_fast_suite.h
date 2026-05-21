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

#include "rom_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosRomFastSuite : public KratosCoreFastSuite {
public:
  KratosRomFastSuite();

private:
  KratosRomApplication::Pointer mpRomApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
