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

#include "trilinos_application.h"
#include "mpi/testing/mpi_testing.h"

namespace Kratos::Testing {

class KratosTrilinosFastSuite : public KratosMPICoreFastSuite {
public:
  KratosTrilinosFastSuite();

private:
  KratosTrilinosApplication::Pointer mpTrilinosApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

class KratosTrilinosApplicationMPITestSuite : public KratosTrilinosFastSuite {};

} // namespace Kratos::Testing
