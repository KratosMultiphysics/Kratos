// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

namespace Kratos
{
    class KratosGeoParser
    {
    public:
        KratosGeoParser() = delete;
        ~KratosGeoParser() {};
        static void parseMaterial(Model& model, std::string filepath);
        static void parseMesh(ModelPart& model_part, std::string filepath);
        static Parameters openProjectParamsFile(std::string filepath);
        static std::vector<std::shared_ptr<Process>> parseProcess(ModelPart& model_part, Parameters projFile);
    };
}