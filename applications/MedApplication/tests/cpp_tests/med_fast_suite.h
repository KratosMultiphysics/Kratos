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

#include "med_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosMedFastSuite : public KratosCoreFastSuite {
public:
  KratosMedFastSuite();

private:
  KratosMeshingApplication::Pointer mpMedApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
