//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "output_quadrature_domain_process.h"

namespace Kratos
{

OutputQuadratureDomainProcess::OutputQuadratureDomainProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

void OutputQuadratureDomainProcess::ExecuteBeforeSolutionLoop()
{
    bool OutputGeometryElements = mThisParameters["output_geometry_elements"].GetBool();
    bool OutputGeometryConditions = mThisParameters["output_geometry_conditions"].GetBool();

    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);
    std::string output_file_name = mThisParameters["output_file_name"].GetString();

    std::string contents = "{ \"geometry_integration_points\":[";

    if (OutputGeometryElements) {
        for (auto element : r_model_part.Elements()) {
            auto integration_point = element.GetGeometry().IntegrationPoints()[0];
            contents += '[' + std::to_string(element.Id()) + ',' + std::to_string(element.GetGeometry().GetGeometryParent(0).Id()) + ",[";
            contents += std::to_string(integration_point[0]) + ',' + std::to_string(integration_point[1]) + "]],";
        }
    }
    if (OutputGeometryConditions) {
        for (auto condition : r_model_part.Conditions()) {
            auto integration_point = condition.GetGeometry().IntegrationPoints()[0];
            contents += '[' + std::to_string(condition.Id()) + ',' + std::to_string(condition.GetGeometry().GetGeometryParent(0).Id()) + ",[";
            contents += std::to_string(integration_point[0]) + ',' + std::to_string(integration_point[1]) + "]],";
        }
    }
    /// cut off last ","
    contents.pop_back();
    contents += "]}";

    std::ofstream file;
    file.open(output_file_name);
    file << contents;
    file.close();
}

const Parameters OutputQuadratureDomainProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "output_file_name"           : "",
        "model_part_name"            : "",
        "output_geometry_elements"   : true,
        "output_geometry_conditions" : false
    })" );
    return default_parameters;
}

} // namespace Kratos
