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

// Project includes
#include "testing/testing.h"
#include "constitutive_laws_application.h"
#include "structural_mechanics_application.h"

namespace Kratos::Testing 
{

class KratosConstitutiveLawsFastSuite : public KratosCoreFastSuite 
{
  public:
    KratosConstitutiveLawsFastSuite();

  private:
    KratosConstitutiveLawsApplication::Pointer mpConstitutiveLawsApp;
};

class KratosConstitutiveLawsWithStructuralElementsSuite : public KratosCoreFastSuite 
{
  public:
    KratosConstitutiveLawsWithStructuralElementsSuite();

  private:
    KratosStructuralMechanicsApplication::Pointer mpStructuralMechanicsApp;
    KratosConstitutiveLawsApplication::Pointer mpConstitutiveLawsApp;
};

} // namespace Kratos::Testing
