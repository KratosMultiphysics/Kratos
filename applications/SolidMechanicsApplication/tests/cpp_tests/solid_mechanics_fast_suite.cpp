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

// External includes
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// Project includes
#include "solid_mechanics_fast_suite.h"

namespace Kratos::Testing{
KratosSolidMechanicsFastSuite::KratosSolidMechanicsFastSuite() 
  : KratosCoreFastSuite() {
  mpSolidMechanicsApp = std::make_shared<KratosSolidMechanicsApplication>();
   this->ImportApplicationIntoKernel(mpSolidMechanicsApp);
}

} // namespace Kratos::Testing


