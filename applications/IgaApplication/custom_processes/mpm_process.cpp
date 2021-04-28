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


// External includes

// Project includes
#include "mpm_process.h"

namespace Kratos
{

MpmProcess::MpmProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

void MpmProcess::ExecuteInitializeSolutionStep()
{
    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);

    for (auto& element : r_model_part.Elements()) {
        auto& r_geometry = element.GetGeometry();

        /// Update location
        array_1d<double, 3> mpm_location;
        r_geometry.GlobalCoordinates(mpm_location, 0);
        array_1d<double, 3> local_coordinates_2 = element.GetGeometry().IntegrationPoints()[0];
        element.SetValue(MPM_LOCATION, mpm_location);
    }
    ResetNodalVariables(r_model_part);

    for (auto& element : r_model_part.Elements()) {
        const auto& mpm_location = element.GetValue(MPM_LOCATION);
        auto& r_parent_geometry = element.GetGeometry().GetGeometryParent(0);
        array_1d<double, 3> local_coordinates = element.GetGeometry().IntegrationPoints()[0];
        array_1d<double, 3> global_coordinates = ZeroVector(3);
        int check = r_parent_geometry.ProjectionPoint(mpm_location, global_coordinates, local_coordinates);
        if (check == 0) {
            element.Set(TO_ERASE);
            //r_model_part.RemoveElementFromAllLevels(element.Id());
            KRATOS_WATCH("Element was set to be removed")
            KRATOS_WATCH(mpm_location)
        }
        else {
            IntegrationPoint<3> integration_point(local_coordinates, element.GetValue(MPM_WEIGHT));
            std::vector<IntegrationPoint<3>> integration_points_array(1);
            PointerVector<Geometry<Node<3>>> geometries(1);
            integration_points_array[0] = integration_point;
            r_parent_geometry.CreateQuadraturePointGeometries(geometries, 3, integration_points_array);
            auto& geom = element.GetGeometry();
            geom = geometries[0];

            element.SetGeometry(geometries(0));
        }
    }

}

void MpmProcess::ExecuteFinalizeSolutionStep()
{
}

void MpmProcess::ResetNodalVariables(ModelPart& rModelPart) {
    VariableUtils().UpdateCurrentToInitialConfiguration(rModelPart.Nodes());

    #pragma omp parallel for
    for (int iter = 0; iter < static_cast<int>(rModelPart.Nodes().size()); ++iter)
    {
        auto i = rModelPart.NodesBegin() + iter;
        array_1d<double, 3 >& r_nodal_displacement = i->FastGetSolutionStepValue(DISPLACEMENT);

        r_nodal_displacement.clear();
    }
}

const Parameters MpmProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name" : ""
    })" );
    return default_parameters;
}

} // namespace Kratos
