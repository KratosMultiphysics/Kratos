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
#include "contact_structural_mechanics_application.h"

namespace Kratos::Testing 
{

class KratosContactStructuralMechanicsFastSuite : public KratosCoreFastSuite 
{
  public:
    KratosContactStructuralMechanicsFastSuite();

  private:
    KratosContactStructuralMechanicsApplication::Pointer mpContactStructuralMechanicsApp;
};

} // namespace Kratos::Testing
