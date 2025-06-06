//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#pragma once

#include "testing/testing.h"

namespace Kratos
{
class KratosIgaApplication;
class KratosStructuralMechanicsApplication;
} // namespace Kratos

namespace Kratos::Testing
{

class KratosIgaSMFastSuite : public KratosCoreFastSuite
{
public:
    KratosIgaSMFastSuite();

private:
    std::shared_ptr<KratosIgaApplication>  mpIgaApp;
    std::shared_ptr<KratosStructuralMechanicsApplication> mpSMApp;
};

} // namespace Kratos::Testing
