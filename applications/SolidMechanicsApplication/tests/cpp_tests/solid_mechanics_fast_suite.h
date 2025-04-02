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
#include "solid_mechanics_application.h"

namespace Kratos::Testing 
{

class SolidMechanicsApplicationFastSuite : public KratosCoreFastSuite {
public:
  SolidMechanicsApplicationFastSuite();

private:
  KratosSolidMechanicsApplication::Pointer mpSolidMechanicsApp;
};

} // namespace Kratos::Testing
