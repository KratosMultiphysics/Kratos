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
    Parameters ThisParameters
    ) : mrModel(rModel)
      , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

void OutputQuadratureDomainProcess::ExecuteBeforeSolutionLoop()
{
    bool OutputGeometryElements = mThisParameters["output_geometry_elements"].GetBool();
    bool OutputGeometryConditions = mThisParameters["output_geometry_conditions"].GetBool();
    bool OutputBrepElements = mThisParameters["output_brep_elements"].GetBool();
    bool OutputBrepConditions = mThisParameters["output_brep_conditions"].GetBool();

    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);
    std::string output_file_name = mThisParameters["output_file_name"].GetString();

    std::ofstream file;
    file.open(output_file_name);

    file << "{ \"geometry_integration_points\":["

    for (auto element : r_model_part.Elements()) {
        auto integration_point = element.GetGeometry().IntegtrationPoints()[0];
        file << "[" << element.Id << "," << element.GetGeometryParent(0).Id() << ",[" << integration_point[0] << "," << integration_point[0] << "]],"
    }
    /// cut off last ","
    file.pop_back();
    file << "}";
    file.close();
}

const Parameters OutputQuadratureDomainProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "output_file_name"           : "",
        "model_part_name"            : "",
        "output_geometry_elements"   : true,
        "output_geometry_conditions" : false,
        "output_brep_elements"       : false,
        "output_brep_conditions"     : false
    })" );
    return default_parameters;
}

} // namespace Kratos
