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
//

#pragma once

#include <string>

#include <includes/kratos_parameters.h>


namespace Kratos
{

class Model;

class InputUtilities
{
public:
    static Parameters ProjectParametersFrom(const std::string& rProjectFilePath);
    static void AddMaterialsFrom(const std::string& rMaterialFilePath,
                                 Model&             rModel);
};

}
