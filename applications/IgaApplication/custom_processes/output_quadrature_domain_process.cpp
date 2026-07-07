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
    bool OutputCouplingGeometryConditions = mThisParameters["output_coupling_geometry_conditions"].GetBool();

    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);
    std::string output_file_name = mThisParameters["output_file_name"].GetString();

    std::string contents = "{\n\"geometry_integration_points\":[ \n";

    if (OutputGeometryElements) {
        for (auto element : r_model_part.Elements()) {
            auto integration_point = element.GetGeometry().IntegrationPoints()[0];
            contents += '[' + std::to_string(element.Id()) + ',' + std::to_string(element.GetGeometry().GetGeometryParent(0).Id()) + ",[";
            contents += std::to_string(integration_point[0]) + ',' + std::to_string(integration_point[1]) + "]],\n";
        }
    }
    if (OutputGeometryConditions) {
        for (auto condition : r_model_part.Conditions()) {
            auto integration_point = condition.GetGeometry().IntegrationPoints()[0];
            contents += '[' + std::to_string(condition.Id()) + ',' + std::to_string(condition.GetGeometry().GetGeometryParent(0).Id()) + ",[";
            contents += std::to_string(integration_point[0]) + ',' + std::to_string(integration_point[1]) + "]],\n";
        }
    }
    /// cut off last ","
    contents.pop_back();
    contents.pop_back();
    contents += "\n]";

    if (OutputCouplingGeometryConditions)
    {
        contents += ",\n\"geometry_coupling_integration_points\":[\n";
        for (auto condition : r_model_part.Conditions()) {
            if (condition.GetGeometry().NumberOfGeometryParts() > 1) {
                const auto p_master = condition.GetGeometry().pGetGeometryPart(CouplingGeometry<Node>::Master);
                const auto p_slave = condition.GetGeometry().pGetGeometryPart(CouplingGeometry<Node>::Slave);

                array_1d<double, 3> local_coordinates_on_master_patch = p_master->IntegrationPoints()[0];
                p_master->GetGeometryParent(0).Calculate(PARAMETER_2D_COORDINATES, local_coordinates_on_master_patch);
                array_1d<double, 3> local_coordinates_on_slave_patch = p_slave->IntegrationPoints()[0];
                p_slave->GetGeometryParent(0).Calculate(PARAMETER_2D_COORDINATES, local_coordinates_on_slave_patch);

                IndexType master_patch_id = p_master->GetGeometryParent(0).GetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX).GetGeometryParent(0).Id();
                IndexType slave_patch_id = p_slave->GetGeometryParent(0).GetGeometryPart(Geometry<Node>::BACKGROUND_GEOMETRY_INDEX).GetGeometryParent(0).Id();

                contents += '[' + std::to_string(condition.Id()) + ',' + std::to_string(master_patch_id) + ",[";
                contents += std::to_string(local_coordinates_on_master_patch[0]) + ',' + std::to_string(local_coordinates_on_master_patch[1]) + "],";
                contents += std::to_string(slave_patch_id) + ",[";
                contents += std::to_string(local_coordinates_on_slave_patch[0]) + ',' + std::to_string(local_coordinates_on_slave_patch[1]) + "]],\n";
            }
        }
        /// cut off last ","
        contents.pop_back();
        contents.pop_back();
        contents += "\n]";
    }

    contents += "\n}";
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
        "output_geometry_conditions" : false,
        "output_coupling_geometry_conditions" : false
    })" );
    return default_parameters;
}

} // namespace Kratos
