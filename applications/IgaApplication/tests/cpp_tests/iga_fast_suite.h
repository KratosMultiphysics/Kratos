//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  Kratos default license: kratos/license.txt
//
//  Main authors:   Carlos Roig
//                  Andrea Gorgi

#pragma once

#include "testing/testing.h"

namespace Kratos
{
class KratosIgaApplication;
} // namespace Kratos

namespace Kratos::Testing
{

class KratosIgaFastSuite : public KratosCoreFastSuite
{
public:
    KratosIgaFastSuite();

private:
    std::shared_ptr<KratosIgaApplication>  mpIgaApp;
};

class KratosIgaFast5PSuite : public KratosCoreFastSuite
{
public:
    KratosIgaFast5PSuite();

private:
    std::shared_ptr<KratosIgaApplication>  mpIgaApp;
};

} // namespace Kratos::Testing