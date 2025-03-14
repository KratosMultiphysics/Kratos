//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dagmawi Bekel
//                   Ruben Zorrilla
//
#if USE_TRIANGLE_NONFREE_TPL
// Project includes
#include "includes/define.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "cad_tessellation_modeler.h"

namespace Kratos
{
///@name Operations
///@{

Modeler::Pointer CadTessellationModeler::Create(Model& rModel,
    const Parameters ModelParameters) const
{
    return Kratos::make_shared<CadTessellationModeler>(rModel, ModelParameters);
}

const Parameters CadTessellationModeler::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "cad_model_part_name": "",
        "skin_model_part_name": "",
        "absolute_chordal_error"   : 1e-2,
        "absolute_triangulation_error"  : 1e-2,
        "initial_triangle_area"         : 1,
        "max_triangulation_iteration"   : 10,
        "echo_level": 0
    })");

    return default_parameters;
}

///@}
///@name Stages
///@{

void CadTessellationModeler::SetupModelPart()
{

    const Parameters default_parameters = this->GetDefaultParameters();
    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
        << "Missing \"cad_model_part_name\" in CadTessellationModeler Parameters" << std::endl;
    ModelPart& cad_model_part =
        mpModel->GetModelPart(mParameters["cad_model_part_name"].GetString());

    KRATOS_ERROR_IF_NOT(mParameters.Has("skin_model_part_name"))
        << "Missing \"skin_model_part_name\" in CadTessellationModeler Parameters" << std::endl;
    const std::string skin_model_part_name = mParameters["skin_model_part_name"].GetString();
    ModelPart &skin_model_part = mpModel->HasModelPart(skin_model_part_name)
                                     ? mpModel->GetModelPart(skin_model_part_name)
                                     : mpModel->CreateModelPart(skin_model_part_name);

    KRATOS_ERROR_IF(skin_model_part.NumberOfNodes() != 0)
        << "Skin model part should be empty. Current number of nodes: "
        << skin_model_part.NumberOfNodes() << std::endl;

    IndexType node_id = 0;
    IndexType element_id = 0;
    IndexType vertex_id = 0;
    const auto& r_geometries = cad_model_part.Geometries();
    for (auto it = r_geometries.begin(); it != r_geometries.end(); ++it) {

        IndexType trim_index = 0;
        if (it->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Brep_Surface) {

            const auto& r_aux_geometry = *it;
            const auto r_brep_surface_geom = dynamic_cast<const BrepSurfaceType&>(r_aux_geometry);

            std::vector<array_1d<double, 2>> boundary_loop_uv;
            IndexType point_id = 0;

            while (r_brep_surface_geom.HasGeometryPart(trim_index)) {

                auto p_surface_geometry = cad_model_part.pGetGeometry(r_brep_surface_geom.Id());
                auto p_surface_trim = p_surface_geometry->pGetGeometryPart(trim_index);
                auto p_brep_curve_on_surface =
                    dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_surface_trim);

                auto boundary_segment_tessellation = this->ComputeBoundaryTessellation(*p_brep_curve_on_surface);

                boundary_loop_uv.resize(boundary_loop_uv.size() + boundary_segment_tessellation.size(), ZeroVector(2));

                for (IndexType i = 0; i < boundary_segment_tessellation.size(); ++i) {

                    boundary_loop_uv[point_id][0] = boundary_segment_tessellation[i][0];
                    boundary_loop_uv[point_id][1] = boundary_segment_tessellation[i][1];

                    point_id++;
                }
                trim_index++;
            }

            auto triangulation_uv = this->ComputeSurfaceTriangulation(
                r_brep_surface_geom,
                boundary_loop_uv);

            array_1d<double,3> result = ZeroVector(3);
            array_1d<double,3> local_coordinate = ZeroVector(3);

            for (IndexType tri_i = 0; tri_i < triangulation_uv.size(); ++tri_i) {
                for (IndexType point_i = 0; point_i < 3; ++point_i) {

                    local_coordinate[0] = triangulation_uv[tri_i](point_i, 0);
                    local_coordinate[1] = triangulation_uv[tri_i](point_i, 1);
                    auto point_xyz = r_brep_surface_geom.GlobalCoordinates(result, local_coordinate);

                    skin_model_part.CreateNewNode(
                        node_id++,
                        point_xyz[0],
                        point_xyz[1],
                        point_xyz[2]);
                }
            }

            auto p_properties = Kratos::make_shared<Properties>(0);

            // create elements in skin_model_part
            for (IndexType element_i = 0; element_i < triangulation_uv.size(); ++element_i)
            {
                skin_model_part.CreateNewElement("Element3D3N",
                    element_id++,
                    {{vertex_id, vertex_id + 1, vertex_id + 2}},
                    p_properties);

                vertex_id += 3;
            }
        }
    }
}

///@}
///@name Private Operations
///@{

std::vector<array_1d<double, 2>> CadTessellationModeler::ComputeBoundaryTessellation(const BrepCurveOnSurfaceType& rBoundarySegment)
{
    const double chordal_error = mParameters["absolute_chordal_error"].GetDouble();

    auto tessellation = CurveTessellation<ContainerNodeType>::ComputeTessellation(
        rBoundarySegment,
        rBoundarySegment.PolynomialDegree(0),
        rBoundarySegment.pGetCurveOnSurface()->pGetCurve()->DomainInterval(),
        rBoundarySegment.pGetCurveOnSurface()->pGetCurve()->KnotSpanIntervals(),
        chordal_error);

    std::vector<array_1d<double, 2>> boundary_segment_uv;
    array_1d<double,3> result = ZeroVector(3);
    array_1d<double,3> local_coordinate = ZeroVector(3);
    boundary_segment_uv.resize(tessellation.size() - 1, ZeroVector(2));

    for (IndexType i = 0; i < tessellation.size() - 1; ++i)
    {
        local_coordinate[0] = tessellation[i].first;

        boundary_segment_uv[i][0] = rBoundarySegment.pGetCurveOnSurface()->pGetCurve()->GlobalCoordinates(result, local_coordinate)[0];
        boundary_segment_uv[i][1] = rBoundarySegment.pGetCurveOnSurface()->pGetCurve()->GlobalCoordinates(result, local_coordinate)[1];
    }
    return boundary_segment_uv;
}


std::vector<BoundedMatrix<double,3,3>> CadTessellationModeler::ComputeSurfaceTriangulation(
    const BrepSurfaceType& rSurfaceGeometry,
    const std::vector<array_1d<double, 2>>& rBoundaryLoop)
{
    const double absolute_triangulation_error = mParameters["absolute_triangulation_error"].GetDouble();
    double aux_area = mParameters["initial_triangle_area"].GetDouble();
    const IndexType max_iteration = mParameters["max_triangulation_iteration"].GetInt();

    // Initialize a 1d list with the coordinates of the points (outer and inner polygons)
    SizeType n_boundary_points = rBoundaryLoop.size();
    std::vector<double> boundary_points_coords(2 * n_boundary_points);
    IndexPartition<IndexType>(n_boundary_points).for_each([&](IndexType iPoint){
        auto& r_coords = rBoundaryLoop[iPoint];
        boundary_points_coords[2 * iPoint] = r_coords[0];
        boundary_points_coords[2 * iPoint + 1] = r_coords[1];
    });

    // Initilize the segment list with the number of boundary edges and the start and end node id
    // For closed polygons the number of segments is equal to the number of points
    const SizeType n_segments = rBoundaryLoop.size();
    std::vector<std::array<double,2>> segments(n_segments);
    IndexPartition<IndexType>(n_segments).for_each(std::array<double,2>(), [&](IndexType iSegment, std::array<double,2>& rSegmentTLS){
        rSegmentTLS[0] = iSegment;
        rSegmentTLS[1] = iSegment + 1 != n_segments ? iSegment + 1 : 0;
        segments[iSegment] = rSegmentTLS;
    });

    std::pair<std::vector<IndexType>, std::vector<double>> delaunator_output;
    for (IndexType i = 1; i < max_iteration + 1; ++i)
    {
        // Call the triangle within the DelaunatorUtilities to calculate the current patch tessellation
        delaunator_output = DelaunatorUtilities::ComputeTrianglesConnectivity(
            boundary_points_coords,
            segments,
            aux_area);

        const auto& r_output_connectivites = std::get<0>(delaunator_output);
        const auto& r_output_coordinates_list = std::get<1>(delaunator_output);

        auto gauss_points_exact_xyz = this->InsertGaussPointsExactSurface(
            rSurfaceGeometry,
            r_output_coordinates_list,
            r_output_connectivites);

        auto gauss_points_approx_xyz = this->InsertGaussPointsApproxSurface(
            rSurfaceGeometry,
            r_output_coordinates_list,
            r_output_connectivites);

        const double it_error = this->ComputeDiscretizationError(
            gauss_points_exact_xyz,
            gauss_points_approx_xyz);

        KRATOS_INFO_IF("CadTessellationModeler", mEchoLevel >= 1) << "Iteration " << i << std::endl;
        KRATOS_INFO_IF("CadTessellationModeler", mEchoLevel >= 1) << "Area: " << aux_area << " - error: " << it_error << " - tolerance: " << absolute_triangulation_error << std::endl;

        // Check if triangulation error is reached and break out of the loop
        if (absolute_triangulation_error > it_error) {
            break;
        }

        // If error is above a certain value -> remesh again with the area halved
        aux_area /= 2;
    }

    const auto& r_output_connectivites = std::get<0>(delaunator_output);
    const auto& r_output_coordinates_list = std::get<1>(delaunator_output);
    const SizeType n_triangles = r_output_connectivites.size() % 3 == 0 ? r_output_connectivites.size() / 3 : KRATOS_ERROR << "Error in connectivities vector size." << std::endl;
    std::vector<BoundedMatrix<double,3,3>> triangulation_uv(n_triangles);
    for (IndexType i_triangle = 0; i_triangle < n_triangles; ++i_triangle) {
        auto& r_triangle_uv_coords = triangulation_uv[i_triangle];
        for (IndexType j = 0; j < 3; ++j) {
            const IndexType node_id = r_output_connectivites[3 * i_triangle + j];
            r_triangle_uv_coords(j,0) = r_output_coordinates_list[node_id * 2];
            r_triangle_uv_coords(j,1) = r_output_coordinates_list[node_id * 2 + 1];
            r_triangle_uv_coords(j,2) = 0.0;
        }
    }

    return triangulation_uv;
}

std::vector<BoundedMatrix<double,3,3>> CadTessellationModeler::InsertGaussPointsExactSurface(
    const BrepSurfaceType& rSurfaceGeometry,
    const std::vector<double>& rPointsCoordinates,
    const std::vector<IndexType>& rTriangleConnectivities)
{
    const auto gp_canonical_tri = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3>>::GenerateIntegrationPoints();

    const SizeType n_triangles = rTriangleConnectivities.size() % 3 == 0 ? rTriangleConnectivities.size() / 3 : KRATOS_ERROR << "Error in connectivities vector size." << std::endl;
    std::vector<BoundedMatrix<double,3,3>> gp_xyz(n_triangles);
    typedef std::pair<array_1d<double,3>, array_1d<double,3>> TLSType;
    IndexPartition<IndexType>(n_triangles).for_each(TLSType(), [&](IndexType iTriangle, TLSType& rTLSContainer){
        auto& r_point_xyz = std::get<0>(rTLSContainer);
        auto& r_local_coordinate = std::get<1>(rTLSContainer);
        r_local_coordinate[2] = 0.0; // Note that we work in the 2D plane only
        IndexType aux = 3 * iTriangle;
        auto& r_gp_xyz_iTriangle = gp_xyz[iTriangle];
        for (IndexType gp_i = 0; gp_i < 3; ++gp_i) {
            // Calculate current Gauss point local coordinates
            const auto aux_gp = gp_canonical_tri[gp_i];
            r_local_coordinate[0] = rPointsCoordinates[rTriangleConnectivities[aux + 0] * 2] * (1 - aux_gp[0] - aux_gp[1]) +
                                  rPointsCoordinates[rTriangleConnectivities[aux + 1] * 2] * aux_gp[0] +
                                  rPointsCoordinates[rTriangleConnectivities[aux + 2] * 2] * aux_gp[1];
            r_local_coordinate[1] = rPointsCoordinates[rTriangleConnectivities[aux + 0] * 2 + 1] * (1 - aux_gp[0] - aux_gp[1]) +
                                  rPointsCoordinates[rTriangleConnectivities[aux + 1] * 2 + 1] * aux_gp[0] +
                                  rPointsCoordinates[rTriangleConnectivities[aux + 2] * 2 + 1] * aux_gp[1];

            // Calculate current Gauss point global coordinates on exact surface
            rSurfaceGeometry.GlobalCoordinates(r_point_xyz, r_local_coordinate);
            r_gp_xyz_iTriangle(gp_i, 0) = r_point_xyz[0];
            r_gp_xyz_iTriangle(gp_i, 1) = r_point_xyz[1];
            r_gp_xyz_iTriangle(gp_i, 2) = r_point_xyz[2];
        }
    });

    return gp_xyz;
}

std::vector<BoundedMatrix<double,3,3>> CadTessellationModeler::InsertGaussPointsApproxSurface(
    const BrepSurfaceType& rSurfaceGeometry,
    const std::vector<double>& rPointsCoordinates,
    const std::vector<IndexType>& rTriangleConnectivities)
{
    const auto gp_canonical_tri = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3>>::GenerateIntegrationPoints();

    const SizeType n_triangles = rTriangleConnectivities.size() % 3 == 0 ? rTriangleConnectivities.size() / 3 : KRATOS_ERROR << "Error in connectivities vector size." << std::endl;
    std::vector<BoundedMatrix<double,3,3>> gp_xyz(n_triangles);
    typedef std::tuple<array_1d<double,3>, array_1d<double,3>, array_1d<double,3>, array_1d<double,3>> TLSType;
    IndexPartition<IndexType>(n_triangles).for_each(TLSType(), [&](IndexType iTriangle, TLSType& rTLSContainer){
        IndexType aux = 3 * iTriangle;
        auto& r_point_0 = std::get<0>(rTLSContainer);
        auto& r_point_1 = std::get<1>(rTLSContainer);
        auto& r_point_2 = std::get<2>(rTLSContainer);
        auto& r_local_coordinate = std::get<3>(rTLSContainer);
        r_local_coordinate[2] = 0.0; // Note that we work in the 2D plane only

        // Calculate triangle point 0 global coordinates over the BREP surface
        r_local_coordinate[0] = rPointsCoordinates[rTriangleConnectivities[aux + 0] * 2];
        r_local_coordinate[1] = rPointsCoordinates[rTriangleConnectivities[aux + 0] * 2 + 1];
        rSurfaceGeometry.GlobalCoordinates(r_point_0, r_local_coordinate);

        // Calculate triangle point 1 global coordinates over the BREP surface
        r_local_coordinate[0] = rPointsCoordinates[rTriangleConnectivities[aux + 1] * 2];
        r_local_coordinate[1] = rPointsCoordinates[rTriangleConnectivities[aux + 1] * 2 + 1];
        rSurfaceGeometry.GlobalCoordinates(r_point_1, r_local_coordinate);

        // Calculate triangle point 2 global coordinates over the BREP surface
        r_local_coordinate[0] = rPointsCoordinates[rTriangleConnectivities[aux + 2] * 2];
        r_local_coordinate[1] = rPointsCoordinates[rTriangleConnectivities[aux + 2] * 2 + 1];
        rSurfaceGeometry.GlobalCoordinates(r_point_2, r_local_coordinate);

        // Save the approximate surface Gauss points coordinates
        auto& r_gp_xyz_iTriangle = gp_xyz[iTriangle];
        for (IndexType gp_i = 0; gp_i < 3; ++gp_i) {
            const auto aux_gp = gp_canonical_tri[gp_i];
            for (IndexType j = 0; j < 3; ++j) {
                r_gp_xyz_iTriangle(gp_i, j) = r_point_0[j] * (1 - aux_gp[0] - aux_gp[1]) + r_point_1[j] * aux_gp[0] + r_point_2[j] * aux_gp[1];
            }
        }
    });

    return gp_xyz;
}

double CadTessellationModeler::ComputeDiscretizationError(
    const std::vector<BoundedMatrix<double,3,3>>& rGaussPointsExact,
    const std::vector<BoundedMatrix<double,3,3>>& rGaussPointsApprox)
{
    const SizeType n_triangles = rGaussPointsExact.size();
    double discretization_error = IndexPartition<IndexType>(n_triangles).for_each<MaxReduction<double>>([&](IndexType iTriangle){
        double max_triangle_error = 0.0;
        const auto& r_gauss_pt_exact = rGaussPointsExact[iTriangle];
        const auto& r_gauss_pt_approx = rGaussPointsApprox[iTriangle];
        for (IndexType i_point = 0; i_point < 3; ++i_point) {
            double point_error = 0.0;
            for (IndexType dim = 0; dim < 3; ++dim) {
                point_error += std::pow(r_gauss_pt_exact(i_point,dim)-r_gauss_pt_approx(i_point,dim), 2);
            }
            point_error = std::sqrt(point_error);
            if (point_error > max_triangle_error) {
                max_triangle_error = point_error;
            }
        }
        return max_triangle_error;
    });

    return discretization_error;
}

} // namespace Kratos
#endif