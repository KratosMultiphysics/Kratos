//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes

// Project includes
#include "custom_processes/prepare_integration_on_true_boundary_process.h"
#include "iga_application_variables.h"
#include "integration/integration_point_utilities.h"

namespace Kratos
{

PrepareIntegrationOnTrueBoundaryProcess::PrepareIntegrationOnTrueBoundaryProcess(
    Model& rModel,
    Parameters ThisParameters)
    : mrModel(rModel)
    , mParameters(ThisParameters)
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

const Parameters PrepareIntegrationOnTrueBoundaryProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "analysis_model_part_name" : "",
        "skin_model_part_name" : "skin_model_part",
        "precision_order_on_integration" : 0
    })");
}

void PrepareIntegrationOnTrueBoundaryProcess::Execute()
{
    const std::string analysis_model_part_name = mParameters["analysis_model_part_name"].GetString();
    KRATOS_ERROR_IF(analysis_model_part_name.empty())
        << "PrepareIntegrationOnTrueBoundaryProcess requires a non-empty \"analysis_model_part_name\"." << std::endl;

    ModelPart& r_analysis_model_part = mrModel.GetModelPart(analysis_model_part_name);
    if (!r_analysis_model_part.HasSubModelPart("SBM_Support_inner")) {
        return;
    }

    const int domain_size = r_analysis_model_part.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(domain_size != 2 && domain_size != 3)
        << "PrepareIntegrationOnTrueBoundaryProcess supports only 2D and 3D SBM inner boundaries." << std::endl;

    ModelPart& r_surrogate_inner_model_part = r_analysis_model_part.GetSubModelPart("SBM_Support_inner");
    if (r_surrogate_inner_model_part.NumberOfConditions() == 0) {
        return;
    }

    const std::string skin_model_part_name = mParameters["skin_model_part_name"].GetString();
    KRATOS_ERROR_IF_NOT(mrModel.HasModelPart(skin_model_part_name))
        << "PrepareIntegrationOnTrueBoundaryProcess could not find skin model part \""
        << skin_model_part_name << "\"." << std::endl;

    ModelPart& r_skin_model_part = mrModel.GetModelPart(skin_model_part_name);
    KRATOS_ERROR_IF_NOT(r_skin_model_part.HasSubModelPart("inner"))
        << "PrepareIntegrationOnTrueBoundaryProcess requires \"" << skin_model_part_name
        << ".inner\" to exist." << std::endl;

    ModelPart& r_skin_inner_model_part = r_skin_model_part.GetSubModelPart("inner");
    KRATOS_ERROR_IF(r_skin_inner_model_part.NumberOfConditions() == 0)
        << "PrepareIntegrationOnTrueBoundaryProcess requires a non-empty inner skin boundary."
        << std::endl;

    KRATOS_INFO("PrepareIntegrationOnTrueBoundaryProcess")
        << "Starting true-boundary integration for model part \""
        << analysis_model_part_name << "\" with "
        << r_skin_inner_model_part.NumberOfConditions()
        << " inner skin conditions and "
        << r_surrogate_inner_model_part.NumberOfConditions()
        << " surrogate support conditions." << std::endl;

    this->ResetIntegrationData(r_surrogate_inner_model_part);

    PointVector support_condition_centers;
    support_condition_centers.reserve(r_surrogate_inner_model_part.NumberOfConditions());
    for (auto& r_condition : r_surrogate_inner_model_part.Conditions()) {
        const auto center = r_condition.GetGeometry().Center();
        support_condition_centers.push_back(
            Kratos::make_intrusive<PointType>(r_condition.Id(), center.X(), center.Y(), center.Z()));
    }

    DynamicBins bins(support_condition_centers.begin(), support_condition_centers.end());
    const SizeType number_of_results = support_condition_centers.size();
    PointVector results(number_of_results);
    DistanceVector distances(number_of_results);

    const int precision_order = mParameters["precision_order_on_integration"].GetInt();
    KRATOS_ERROR_IF(precision_order < 0)
        << "\"precision_order_on_integration\" must be non-negative." << std::endl;

    const double search_radius = this->ComputeSearchRadius(r_analysis_model_part);

    SizeType total_assigned_integration_points = 0;
    SizeType touched_surrogate_conditions = 0;

    auto assign_to_nearest_surrogate =
        [&](const array_1d<double, 3>& rTrueBoundaryPoint,
            const double IntegrationWeight,
            const array_1d<double, 3>& rNormal)
    {
        PointType point_to_search(0, rTrueBoundaryPoint[0], rTrueBoundaryPoint[1], rTrueBoundaryPoint[2]);
        const SizeType obtained_results = bins.SearchInRadius(
            point_to_search,
            search_radius,
            results.begin(),
            distances.begin(),
            number_of_results);

        KRATOS_ERROR_IF(obtained_results == 0)
            << "PrepareIntegrationOnTrueBoundaryProcess found no surrogate support condition within radius "
            << search_radius << " for point " << point_to_search << '.' << std::endl;

        SizeType nearest_result_index = 0;
        double minimum_distance = std::numeric_limits<double>::max();
        for (SizeType i = 0; i < obtained_results; ++i) {
            if (distances[i] < minimum_distance) {
                minimum_distance = distances[i];
                nearest_result_index = i;
            }
        }

        Condition& r_surrogate_condition =
            r_surrogate_inner_model_part.GetCondition(results[nearest_result_index]->Id());

        if (r_surrogate_condition.GetValue(INTEGRATION_WEIGHTS).size() == 0) {
            ++touched_surrogate_conditions;
        }

        this->AppendIntegrationData(
            r_surrogate_condition,
            rTrueBoundaryPoint,
            IntegrationWeight,
            rNormal);

        ++total_assigned_integration_points;
    };

    for (auto& r_skin_condition : r_skin_inner_model_part.Conditions()) {
        const auto& r_geometry = r_skin_condition.GetGeometry();
        if (domain_size == 2) {
            const SizeType number_of_integration_points = static_cast<SizeType>(precision_order + 1);
            KRATOS_ERROR_IF(number_of_integration_points == 0)
                << "At least one integration point is required on the true boundary." << std::endl;
            KRATOS_ERROR_IF(number_of_integration_points > IntegrationPointUtilities::s_gauss_legendre.size())
                << "\"precision_order_on_integration\" = " << precision_order
                << " exceeds the available Gauss-Legendre rules." << std::endl;

            KRATOS_ERROR_IF_NOT(r_geometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear &&
                                r_geometry.PointsNumber() == 2)
                << "PrepareIntegrationOnTrueBoundaryProcess expects 2-node inner skin conditions in 2D."
                << " Condition " << r_skin_condition.Id() << " has geometry family "
                << static_cast<int>(r_geometry.GetGeometryFamily()) << " and "
                << r_geometry.PointsNumber() << " nodes." << std::endl;

            const auto& r_integration_rule = IntegrationPointUtilities::s_gauss_legendre[number_of_integration_points - 1];
            const auto& r_start = r_geometry[0].Coordinates();
            array_1d<double, 3> segment = r_geometry[1].Coordinates() - r_start;
            const double segment_length = norm_2(segment);
            KRATOS_ERROR_IF(segment_length < 1e-14)
                << "PrepareIntegrationOnTrueBoundaryProcess found a zero-length inner skin condition: "
                << r_skin_condition.Id() << std::endl;

            array_1d<double, 3> tangent = segment;
            tangent /= segment_length;

            array_1d<double, 3> normal;
            normal[0] = tangent[1];
            normal[1] = -tangent[0];
            normal[2] = 0.0;

            for (const auto& r_gauss_point : r_integration_rule) {
                array_1d<double, 3> true_boundary_point = r_start + segment * r_gauss_point[0];
                const double integration_weight = r_gauss_point[1] * segment_length;
                assign_to_nearest_surrogate(true_boundary_point, integration_weight, normal);
            }
        } else {
            KRATOS_ERROR_IF(precision_order >= static_cast<int>(IntegrationPointUtilities::s_gauss_triangle.size()))
                << "\"precision_order_on_integration\" = " << precision_order
                << " exceeds the available triangle Gauss rules." << std::endl;

            KRATOS_ERROR_IF_NOT(r_geometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle)
                << "PrepareIntegrationOnTrueBoundaryProcess expects triangular inner skin conditions in 3D."
                << " Condition " << r_skin_condition.Id() << " has geometry family "
                << static_cast<int>(r_geometry.GetGeometryFamily()) << '.' << std::endl;

            const auto& r_integration_rule = IntegrationPointUtilities::s_gauss_triangle[precision_order];
            for (const auto& r_gauss_point : r_integration_rule) {
                array_1d<double, 3> local_coordinates = ZeroVector(3);
                local_coordinates[0] = r_gauss_point[0];
                local_coordinates[1] = r_gauss_point[1];

                array_1d<double, 3> true_boundary_point = ZeroVector(3);
                r_geometry.GlobalCoordinates(true_boundary_point, local_coordinates);

                const double integration_weight =
                    r_gauss_point[2] * std::abs(r_geometry.DeterminantOfJacobian(local_coordinates));
                KRATOS_ERROR_IF(integration_weight < 1e-14)
                    << "PrepareIntegrationOnTrueBoundaryProcess found a degenerate inner skin triangle: "
                    << r_skin_condition.Id() << std::endl;

                const array_1d<double, 3> normal = r_geometry.UnitNormal(local_coordinates);
                assign_to_nearest_surrogate(true_boundary_point, integration_weight, normal);
            }
        }
    }

    KRATOS_INFO("PrepareIntegrationOnTrueBoundaryProcess")
        << "True-boundary integration completed: assigned "
        << total_assigned_integration_points << " integration points to "
        << touched_surrogate_conditions << " surrogate support conditions."
        << std::endl;
}

double PrepareIntegrationOnTrueBoundaryProcess::ComputeSearchRadius(
    const ModelPart& rAnalysisModelPart) const
{
    Vector knot_span_sizes;
    if (rAnalysisModelPart.GetParentModelPart().Has(KNOT_SPAN_SIZES)) {
        knot_span_sizes = rAnalysisModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    } else if (rAnalysisModelPart.Has(KNOT_SPAN_SIZES)) {
        knot_span_sizes = rAnalysisModelPart.GetValue(KNOT_SPAN_SIZES);
    } else if (rAnalysisModelPart.GetRootModelPart().Has(KNOT_SPAN_SIZES)) {
        knot_span_sizes = rAnalysisModelPart.GetRootModelPart().GetValue(KNOT_SPAN_SIZES);
    }

    KRATOS_ERROR_IF(knot_span_sizes.size() < 2)
        << "PrepareIntegrationOnTrueBoundaryProcess requires KNOT_SPAN_SIZES with at least two entries."
        << std::endl;

    const int domain_size = rAnalysisModelPart.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(knot_span_sizes.size() < static_cast<SizeType>(domain_size))
        << "PrepareIntegrationOnTrueBoundaryProcess requires KNOT_SPAN_SIZES with at least "
        << domain_size << " entries in " << domain_size << "D." << std::endl;

    double reference_mesh_size = knot_span_sizes[0];
    for (int i = 1; i < domain_size; ++i) {
        reference_mesh_size = std::max(reference_mesh_size, knot_span_sizes[i]);
    }

    return 2.0 * std::sqrt(static_cast<double>(domain_size)) * reference_mesh_size;
}

void PrepareIntegrationOnTrueBoundaryProcess::ResetIntegrationData(ModelPart& rSurrogateModelPart) const
{
    for (auto& r_condition : rSurrogateModelPart.Conditions()) {
        r_condition.SetValue(INTEGRATION_POINTS, Matrix(0, 3));
        r_condition.SetValue(INTEGRATION_WEIGHTS, Vector(0));
        r_condition.SetValue(INTEGRATION_NORMALS, Matrix(0, 3));
    }
}

void PrepareIntegrationOnTrueBoundaryProcess::AppendIntegrationData(
    Condition& rCondition,
    const array_1d<double, 3>& rIntegrationPoint,
    const double IntegrationWeight,
    const array_1d<double, 3>& rNormal) const
{
    Matrix integration_points = rCondition.GetValue(INTEGRATION_POINTS);
    const SizeType old_number_of_points = integration_points.size1();
    integration_points.resize(old_number_of_points + 1, 3, true);
    for (IndexType i = 0; i < 3; ++i) {
        integration_points(old_number_of_points, i) = rIntegrationPoint[i];
    }
    rCondition.SetValue(INTEGRATION_POINTS, integration_points);

    Vector integration_weights = rCondition.GetValue(INTEGRATION_WEIGHTS);
    const SizeType old_number_of_weights = integration_weights.size();
    integration_weights.resize(old_number_of_weights + 1, true);
    integration_weights[old_number_of_weights] = IntegrationWeight;
    rCondition.SetValue(INTEGRATION_WEIGHTS, integration_weights);

    Matrix integration_normals = rCondition.GetValue(INTEGRATION_NORMALS);
    const SizeType old_number_of_normals = integration_normals.size1();
    integration_normals.resize(old_number_of_normals + 1, 3, true);
    for (IndexType i = 0; i < 3; ++i) {
        integration_normals(old_number_of_normals, i) = rNormal[i];
    }
    rCondition.SetValue(INTEGRATION_NORMALS, integration_normals);
}

}  // namespace Kratos
