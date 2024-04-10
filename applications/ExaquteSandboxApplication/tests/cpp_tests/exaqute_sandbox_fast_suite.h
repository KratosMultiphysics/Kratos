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

#include "exaqute_sandbox_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosExaquteSanboxFastSuite : public KratosCoreFastSuite {
public:
  KratosExaquteSanboxFastSuite();

private:
  KratosExaquteSanboxApplication::Pointer mpExaquteSanboxApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
