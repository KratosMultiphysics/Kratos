//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#include "iga_structural_mechanics_fast_suite.h"
#include "iga_application.h"
#include "structural_mechanics_application.h"

namespace Kratos::Testing
{

KratosIgaSMFastSuite::KratosIgaSMFastSuite() : KratosCoreFastSuite()
{
    mpSMApp = std::make_shared<KratosStructuralMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpSMApp);
    mpIgaApp = std::make_shared<KratosIgaApplication>();
    this->ImportApplicationIntoKernel(mpIgaApp);
}

} // namespace Kratos::Testing
