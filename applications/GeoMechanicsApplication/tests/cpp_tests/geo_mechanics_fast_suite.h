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

#include "geo_mechanics_application.h"
#include "linear_solvers_application.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

 class KratosGeoMechanicsFastSuite : public KratosCoreFastSuite
 {
 public:
     KratosGeoMechanicsFastSuite();

 private:
     KratosGeoMechanicsApplication::Pointer mpGeoApp;
     KratosLinearSolversApplication::Pointer mpLinearSolversApp;
 };

 class KratosGeoMechanicsIntegrationSuite : public KratosCoreFastSuite
 {
 public:
     KratosGeoMechanicsIntegrationSuite();

 private:
     KratosGeoMechanicsApplication::Pointer mpGeoApp;
     KratosLinearSolversApplication::Pointer mpLinearSolversApp;
 };

} // namespace Kratos::Testing
