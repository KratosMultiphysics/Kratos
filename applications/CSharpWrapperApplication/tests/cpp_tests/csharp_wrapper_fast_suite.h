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

#include "csharp_wrapper_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosCSharpWrapperFastSuite : public KratosCoreFastSuite {
public:
  KratosCSharpWrapperFastSuite();

private:
  KratosCSharpWrapperApplication::Pointer mpCSharpWrapperlApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
