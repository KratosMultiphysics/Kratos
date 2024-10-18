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

// Project includes
#include "constitutive_laws_fast_suite.h"

namespace Kratos::Testing 
{

KratosConstitutiveLawsFastSuite::KratosConstitutiveLawsFastSuite()
    : KratosCoreFastSuite() 
{
    mpConstitutiveLawsApp = std::make_shared<KratosConstitutiveLawsApplication>();
    this->ImportApplicationIntoKernel(mpConstitutiveLawsApp);
}

KratosConstitutiveLawsWithStructuralElementsSuite::KratosConstitutiveLawsWithStructuralElementsSuite()
    : KratosCoreFastSuite() 
{
    mpStructuralMechanicsApp = std::make_shared<KratosStructuralMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpStructuralMechanicsApp);

    mpConstitutiveLawsApp = std::make_shared<KratosConstitutiveLawsApplication>();
    this->ImportApplicationIntoKernel(mpConstitutiveLawsApp);
}

} // namespace Kratos::Testing