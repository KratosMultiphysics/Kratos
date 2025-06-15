//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Carlos Roig
//                  Andrea Gorgi

#include "iga_fast_suite.h"
#include "iga_application.h"

namespace Kratos::Testing
{

KratosIgaFastSuite::KratosIgaFastSuite() : KratosCoreFastSuite()
{
    mpIgaApp = std::make_shared<KratosIgaApplication>();
    this->ImportApplicationIntoKernel(mpIgaApp);
}

KratosIgaFast5PSuite::KratosIgaFast5PSuite() : KratosCoreFastSuite() {
    mpIgaApp = std::make_shared<KratosIgaApplication>();
    this->ImportApplicationIntoKernel(mpIgaApp);
}

} // namespace Kratos::Testing