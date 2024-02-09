// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

namespace Kratos
{

class KratosGeoSettlement;

class CustomWorkflowFactory
{
public:
    static KratosGeoSettlement* CreateKratosGeoSettlement();
};

} // namespace Kratos
