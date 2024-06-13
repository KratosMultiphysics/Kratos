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

#include "includes/gid_io.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

class Model;
class ModelPart;

class GeoOutputWriter
{
public:
    GeoOutputWriter(const Parameters& rGidOutputSettings, const std::string& rWorkingDirectory, ModelPart& rModelPart);

    void WriteGiDOutput(ModelPart& rModelPart, Parameters Settings, bool WriteHydraulicHeadToNodes = true);

    void FinalizeResults();

private:
    void WriteNodalOutput(const std::vector<std::string>& rOutputItemNames, const ModelPart& rModelPart);
    void WriteIntegrationPointOutput(const std::vector<std::string>& rOutputItemNames, const ModelPart& rModelPart);
    void CalculateNodalHydraulicHead(ModelPart& rModelPart);

    static GidIO<> MakeGidIO(const std::string& rWorkingDirectory, const Parameters& rGidOutputSettings);

    GidIO<> mGidIO;
};

} // namespace Kratos
