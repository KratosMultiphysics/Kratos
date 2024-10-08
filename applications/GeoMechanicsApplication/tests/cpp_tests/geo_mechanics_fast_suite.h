// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Richard Faasse
//

#pragma once

#include "testing/testing.h"

namespace Kratos
{
class KratosGeoMechanicsApplication;
class KratosLinearSolversApplication;
} // namespace Kratos

namespace Kratos::Testing
{

class KratosGeoMechanicsFastSuite : public KratosCoreFastSuite
{
public:
    KratosGeoMechanicsFastSuite();

private:
    std::shared_ptr<KratosGeoMechanicsApplication>  mpGeoApp;
    std::shared_ptr<KratosLinearSolversApplication> mpLinearSolversApp;
};

class KratosGeoMechanicsFastSuiteWithoutKernel : public KratosCoreFastSuiteWithoutKernel
{
public:
    KratosGeoMechanicsFastSuiteWithoutKernel();
};

class KratosGeoMechanicsIntegrationSuite : public KratosCoreFastSuite
{
public:
    KratosGeoMechanicsIntegrationSuite();

private:
    std::shared_ptr<KratosGeoMechanicsApplication>  mpGeoApp;
    std::shared_ptr<KratosLinearSolversApplication> mpLinearSolversApp;
};

} // namespace Kratos::Testing
