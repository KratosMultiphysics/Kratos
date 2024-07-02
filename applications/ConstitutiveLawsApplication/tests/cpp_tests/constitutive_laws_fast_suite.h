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

#include "constitutive_laws_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosConstitutiveLawsFastSuite : public KratosCoreFastSuite {
public:
  KratosConstitutiveLawsFastSuite();

private:
  KratosConstitutiveLawsApplication::Pointer mpConstitutiveLawsApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
