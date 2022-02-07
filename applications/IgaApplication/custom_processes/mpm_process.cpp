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
#include "utilities/quadrature_points_utility.h"

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
        double mpm_weight = r_geometry.IntegrationPoints()[0].Weight()
            * r_geometry.DeterminantOfJacobian(0);
        element.SetValue(MPM_LOCATION, mpm_location);
        element.SetValue(MPM_WEIGHT, mpm_weight);
        //TODO Edge integration.
        
        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry)
        {
            auto normal = r_geometry.Normal(0);
            KRATOS_WATCH(normal)
            r_geometry.Calculate(LOCAL_TANGENT, normal);
        }
    }
    for (auto& element : r_model_part.Conditions()) {
        auto& r_geometry = element.GetGeometry();

        /// Update location
        array_1d<double, 3> mpm_location;
        r_geometry.GlobalCoordinates(mpm_location, 0);
        double mpm_weight = r_geometry.IntegrationPoints()[0].Weight()
            * r_geometry.DeterminantOfJacobian(0);
        element.SetValue(MPM_LOCATION, mpm_location);
        element.SetValue(MPM_WEIGHT, mpm_weight);
        //TODO Edge integration.

        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry)
        {
            auto normal = r_geometry.Normal(0);
            KRATOS_WATCH(normal);
            r_geometry.Calculate(LOCAL_TANGENT, normal);
        }
    }
    ResetNodalVariables(r_model_part);

    auto p_background_geometry = r_model_part.pGetGeometry("background_geometry");
    NurbsSurfaceGeometryPointerType p_nurbs_surface = dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_background_geometry);

    for (auto& element : r_model_part.Elements()) {
        const auto& mpm_location = element.GetValue(MPM_LOCATION);
        const double mpm_weight = element.GetValue(MPM_WEIGHT);
        //auto& r_parent_geometry = element.GetGeometry().GetGeometryParent(0);
        array_1d<double, 3> local_coordinates = element.GetGeometry().IntegrationPoints()[0];
        //array_1d<double, 3> global_coordinates = ZeroVector(3);
        int check = p_nurbs_surface->ProjectionPointGlobalToLocalSpace(
            mpm_location, local_coordinates, 1e-6);
        if (check == 0) {
            element.Set(TO_ERASE);
            //r_model_part.RemoveElementFromAllLevels(element.Id());
            KRATOS_WATCH("Element was set to be removed")
            KRATOS_WATCH(mpm_location)
            KRATOS_WATCH(local_coordinates)
        }
        else {
            auto& r_geometry = element.GetGeometry();
            array_1d<double, 3> r_global_location = r_geometry.Center();
            double weight = r_geometry.IntegrationPoints()[0].Weight() * r_geometry.DeterminantOfJacobian(0);
            IntegrationPoint<3> integration_point(local_coordinates, weight);
            p_nurbs_surface->UpdateQuadraturePointGeometries(
                element.pGetGeometry(), 2, integration_point);
        }
    }
    for (auto& element : r_model_part.Conditions()) {
        const auto& mpm_location = element.GetValue(MPM_LOCATION);
        const double mpm_weight = element.GetValue(MPM_WEIGHT);
        //auto& r_parent_geometry = element.GetGeometry().GetGeometryParent(0);
        array_1d<double, 3> local_coordinates = element.GetGeometry().IntegrationPoints()[0];
        //array_1d<double, 3> global_coordinates = ZeroVector(3);
        int check = p_nurbs_surface->ProjectionPointGlobalToLocalSpace(
            mpm_location, local_coordinates, 1e-6);
        if (check == 0) {
            element.Set(TO_ERASE);
            //r_model_part.RemoveElementFromAllLevels(element.Id());
            KRATOS_WATCH("Element was set to be removed")
                KRATOS_WATCH(mpm_location)
                KRATOS_WATCH(local_coordinates)
        }
        else {
            auto& r_geometry = element.GetGeometry();
            array_1d<double, 3> r_global_location = r_geometry.Center();
            double weight = r_geometry.IntegrationPoints()[0].Weight() * r_geometry.DeterminantOfJacobian(0);
            IntegrationPoint<3> integration_point(local_coordinates, weight);
            p_nurbs_surface->UpdateQuadraturePointGeometries(
                element.pGetGeometry(), 2, integration_point);
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
