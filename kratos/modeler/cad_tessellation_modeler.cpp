//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Project includes
#include "cad_tessellation_modeler.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void CadTessellationModeler::SetupModelPart()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part_name\" in CadTessellationModeler Parameters" << std::endl;
        ModelPart& cad_model_part =
            mpModel->GetModelPart(mParameters["cad_model_part_name"].GetString());

        KRATOS_ERROR_IF_NOT(mParameters.Has("skin_model_part_name"))
            << "Missing \"skin_model_part_name\" in CadTessellationModeler Parameters" << std::endl;
        const std::string skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        ModelPart& skin_model_part = mpModel->HasModelPart(skin_model_part_name)
            ? mpModel->GetModelPart(skin_model_part_name)
            :mpModel->CreateModelPart(skin_model_part_name);

        IndexType node_id = 0;
        IndexType element_id = 0;
        IndexType vertex_id = 0;
        const auto& r_geometries = cad_model_part.Geometries();
        for (auto it = r_geometries.begin(); it != r_geometries.end(); ++it) {

            IndexType trim_index = 0;
            if (it->GetGeometryType() == GeometryData::Kratos_Brep_Surface) {

                const auto& r_aux_geometry = *it;
                const auto r_brep_surface_geom = dynamic_cast<const BrepSurface<PointerVector<Node<3>>, PointerVector<Point>>&>(r_aux_geometry);

                std::vector<array_1d<double, 2>> boundary_loop_uv;
                IndexType point_id = 0;

                while (r_brep_surface_geom.HasGeometryPart(trim_index)) {

                    auto p_surface_geometry = cad_model_part.pGetGeometry(r_brep_surface_geom.Id());
                    auto p_surface_trim = p_surface_geometry->pGetGeometryPart(trim_index);
                    auto p_brep_curve_on_surface =
                        dynamic_pointer_cast<BrepCurveOnSurface<PointerVector<Node<3>>, PointerVector<Point>>>(p_surface_trim);

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

                Properties::Pointer p_properties(new Properties(0));

                // create elements in skin_model_part
                for (IndexType element_i = 0; element_i < triangulation_uv.size(); ++element_i)
                {
                    skin_model_part.CreateNewElement("Element3D3N",
                        element_id++, {{vertex_id, vertex_id + 1, vertex_id + 2}}, p_properties);
                    vertex_id += 3;
                }
            }
        }
    }

    /**
     * @brief This method returns the tessellation of the boundary curves
     */
    std::vector<array_1d<double, 2>> CadTessellationModeler::ComputeBoundaryTessellation(
        const BrepCurveOnSurface<Kratos::Element::NodesArrayType, Kratos::PointerVector<Kratos::Point>>& rBoundaryCurveOnSurface
    )
    {
        auto tessellation = NurbsCurveTessellation<2, PointerVector<Node<3>>>::ComputeTessellation(
            rBoundaryCurveOnSurface,
            rBoundaryCurveOnSurface.PolynomialDegree(0),
            rBoundaryCurveOnSurface.DomainInterval(),
            rBoundaryCurveOnSurface.KnotSpanIntervals(),
            1e-2);

        std::vector<array_1d<double, 2>> boundary_segment_uv;
        array_1d<double,3> result = ZeroVector(3);
        array_1d<double,3> local_coordinate = ZeroVector(3);
        boundary_segment_uv.resize(tessellation.size() - 1, ZeroVector(2));

        for (IndexType i = 0; i < tessellation.size() - 1; ++i)
        {
            local_coordinate[0] = tessellation.at(i).first;

            boundary_segment_uv[i][0] = rBoundaryCurveOnSurface.pGetCurveOnSurface()->pGetCurve()->GlobalCoordinates(result, local_coordinate)[0];
            boundary_segment_uv[i][1] = rBoundaryCurveOnSurface.pGetCurveOnSurface()->pGetCurve()->GlobalCoordinates(result, local_coordinate)[1];
        }
        return boundary_segment_uv;
    }

    /**
     * @brief This method returns the triangulation of a NURBS surface
     */
    std::vector<Matrix> CadTessellationModeler::ComputeSurfaceTriangulation(
        const BrepSurface<Kratos::Element::NodesArrayType, Kratos::PointerVector<Kratos::Point>>& rSurfaceGeometry,
        const std::vector<array_1d<double, 2>>& rBoundaryLoopUV
    )
    {
        // initiliazing the i/o containers
        struct triangulateio in_data;
        struct triangulateio out_data;
        struct triangulateio vor_out_data;

        InitTriangulationDataStructure(in_data);

        // Initialize the pointlist (1d list) with the number of points and the coordinates
        // of the points (outer and inner polygons)
        in_data.numberofpoints = rBoundaryLoopUV.size();
        in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));
        in_data.pointmarkerlist = (int*) malloc(in_data.numberofpoints * sizeof(int));

        IndexType point_id = 0;
        IndexType point_marker_id = 0;
        IndexType point_marker = 0;

        for (IndexType node_i = 0; node_i < rBoundaryLoopUV.size(); ++node_i)
        {
            for (IndexType coords_i = 0; coords_i < 2; ++coords_i)
            {
                in_data.pointlist[point_id++] = rBoundaryLoopUV[node_i][coords_i];
            }
            in_data.pointmarkerlist[point_marker_id++] = point_marker;
        }

        // Initilize the segment list with the number of boundary edges and the start and end node id
        // For closed polygons the number of segments is equal to the number of points
        in_data.numberofsegments = rBoundaryLoopUV.size();
        in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
        in_data.segmentmarkerlist = (int*) malloc(in_data.numberofsegments * sizeof(int));
        IndexType node_id = 0;
        IndexType seg_marker = 0;
        IndexType start_node_id = 0;
        IndexType end_node_id = rBoundaryLoopUV.size();

        for (IndexType seg_i = start_node_id * 2 ; seg_i < end_node_id * 2; ++seg_i)
        {

            in_data.segmentlist[seg_i] = node_id;

            if (node_id == end_node_id)
            {
                in_data.segmentlist[seg_i] = start_node_id;
            }
            if (seg_i % 2 == 0)
            {
                in_data.segmentmarkerlist[seg_i/2] = seg_marker;
                node_id++;
            }
        }

        in_data.numberofholes = 0;
        in_data.holelist = (REAL*) malloc(in_data.numberofholes * 2 * sizeof(REAL));

        auto aux_area = 5.0;
        std::vector<Matrix> triangulation_uv;

        for (IndexType i = 1; i < 2; ++i)
        {
            InitTriangulationDataStructure(out_data);
            InitTriangulationDataStructure(vor_out_data);

            std::string str = "Qqpza" + std::to_string(aux_area);
            char *meshing_options = new char[str.length() + 1];
            strcpy(meshing_options, str.c_str());
            triangulate(meshing_options, &in_data, &out_data, &vor_out_data);

            triangulation_uv.resize(out_data.numberoftriangles, ZeroMatrix(3,2));
            IndexType tri_id = 0;
            for (IndexType i = 0; i < out_data.numberoftriangles; ++i)
            {
                for (IndexType j = 0; j < 3; ++j)
                {
                    triangulation_uv[i](j,0) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2];
                    triangulation_uv[i](j,1) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2 + 1];
                }
                tri_id += 3;
            }
            auto gauss_points_exact_xyz = this->InsertGaussPointsExactSurface(
                rSurfaceGeometry,
                triangulation_uv);

            auto gauss_points_approx_xyz = this->InsertGaussPointsApproxSurface(
                rSurfaceGeometry,
                triangulation_uv);

            auto discretization_error = this->ComputeDiscretizationError(
                gauss_points_exact_xyz,
                gauss_points_approx_xyz);

            auto max_error = *std::max_element(
                std::begin(discretization_error), std::end(discretization_error));
            auto tolerance = false;
            auto max_triangulation_error = 1e-1;

            for (IndexType j = 0; j < discretization_error.size(); ++j)
            {
                if (discretization_error[j] > max_triangulation_error)
                {
                    tolerance = true;
                    break;
                }
            }

            KRATOS_INFO_IF("EMBEDDED_IGA", mEchoLevel >= 0) << "Iteration " << i << std::endl;
            KRATOS_INFO_IF("EMBEDDED_IGA", mEchoLevel >= 0) << "Area: " << aux_area << " - max_error: " << max_error << std::endl;


            // Triangle copies the pointer for the holelist from the in_data to the out_data
            // In order to avoid the freeing of memory twice, which leads to an error the points will
            // be deleted from the out_data and in_data cleans the memory fully.
            out_data.holelist = nullptr;

            CleanTriangulationDataStructure(out_data);
            CleanTriangulationDataStructure(vor_out_data);

            // if error is above a certain value -> remesh again with the area halved
            aux_area /= 2;
        }

        CleanTriangulationDataStructure(in_data);

        return triangulation_uv;
    }

    std::vector<Matrix> CadTessellationModeler::InsertGaussPointsExactSurface(
        const BrepSurface<Kratos::Element::NodesArrayType, Kratos::PointerVector<Kratos::Point>>& rSurfaceGeometry,
        const std::vector<Matrix>& rTriangulation_uv
    )
    {
        /**
         * This function inserts Gauss-Legendre Points into the triangulation of the surface in the parametric domain
         * and subsequently projects these points onto the NURBS surface.
        */

        const auto gp_canonical_tri =
            Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

        array_1d<double,3> result = ZeroVector(3);
        array_1d<double,3> local_coordinate = ZeroVector(3);

        std::vector<Matrix> gp_xyz;
        gp_xyz.resize(rTriangulation_uv.size(), ZeroMatrix(3, 3));


        for (IndexType tri_i = 0; tri_i < rTriangulation_uv.size(); ++tri_i)
        {
            for (IndexType gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i)
            {
                local_coordinate[0] = rTriangulation_uv[tri_i](0,0) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) +
                                      rTriangulation_uv[tri_i](1,0) * gp_canonical_tri[gp_i][0] +
                                      rTriangulation_uv[tri_i](2,0) * gp_canonical_tri[gp_i][1];

                local_coordinate[1] = rTriangulation_uv[tri_i](0,1) * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) +
                                      rTriangulation_uv[tri_i](1,1) * gp_canonical_tri[gp_i][0] +
                                      rTriangulation_uv[tri_i](2,1) * gp_canonical_tri[gp_i][1];

                auto point_xyz = rSurfaceGeometry.GlobalCoordinates(result, local_coordinate);

                gp_xyz[tri_i](gp_i, 0) = point_xyz[0];
                gp_xyz[tri_i](gp_i, 1) = point_xyz[1];
                gp_xyz[tri_i](gp_i, 2) = point_xyz[2];
            }
        }
        return gp_xyz;
    }

    std::vector<Matrix> CadTessellationModeler::InsertGaussPointsApproxSurface(
        const BrepSurface<Kratos::Element::NodesArrayType, Kratos::PointerVector<Kratos::Point>>& rSurfaceGeometry,
        const std::vector<Matrix>& rTriangulation_uv
    )
    {
        /**
         * This function first maps the triangulation into the cartesian space and consequently inserts
         * Gauss-Legendre points into the triangles approximating the exact Surface.
         * These points can be used to measure the distance between the approximated and exact surface
        */

        const auto gp_canonical_tri =
            Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();

        array_1d<double,3> result = ZeroVector(3);
        array_1d<double,3> local_coordinate_0 = ZeroVector(3);
        array_1d<double,3> local_coordinate_1 = ZeroVector(3);
        array_1d<double,3> local_coordinate_2 = ZeroVector(3);

        std::vector<Matrix> gp_xyz;
        gp_xyz.resize(rTriangulation_uv.size(), ZeroMatrix(3, 3));

        for (IndexType tri_i = 0; tri_i < rTriangulation_uv.size(); ++tri_i) {
            for (IndexType gp_i = 0; gp_i < gp_canonical_tri.size(); ++gp_i) {

                local_coordinate_0[0] = rTriangulation_uv[tri_i](0, 0);
                local_coordinate_0[1] = rTriangulation_uv[tri_i](0, 1);
                auto point_0_xyz = rSurfaceGeometry.GlobalCoordinates(result, local_coordinate_0);

                local_coordinate_1[0] = rTriangulation_uv[tri_i](1, 0);
                local_coordinate_1[1] = rTriangulation_uv[tri_i](1, 1);
                auto point_1_xyz = rSurfaceGeometry.GlobalCoordinates(result, local_coordinate_1);

                local_coordinate_2[0] = rTriangulation_uv[tri_i](2, 0);
                local_coordinate_2[1] = rTriangulation_uv[tri_i](2, 1);
                auto point_2_xyz = rSurfaceGeometry.GlobalCoordinates(result, local_coordinate_2);

                gp_xyz[tri_i](gp_i, 0) = point_0_xyz[0] * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) +
                                         point_1_xyz[0] * gp_canonical_tri[gp_i][0] +
                                         point_2_xyz[0] * gp_canonical_tri[gp_i][1];

                gp_xyz[tri_i](gp_i, 1) = point_0_xyz[1] * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) +
                                         point_1_xyz[1] * gp_canonical_tri[gp_i][0] +
                                         point_2_xyz[1] * gp_canonical_tri[gp_i][1];

                gp_xyz[tri_i](gp_i, 2) = point_0_xyz[2] * (1 - gp_canonical_tri[gp_i][0] - gp_canonical_tri[gp_i][1]) +
                                         point_1_xyz[2] * gp_canonical_tri[gp_i][0] +
                                         point_2_xyz[2] * gp_canonical_tri[gp_i][1];

            }
        }
        return gp_xyz;
    }


    Vector CadTessellationModeler::ComputeDiscretizationError(
        const std::vector<Matrix>& rGaussPointsExact,
        const std::vector<Matrix>& rGaussPointsApprox
    )
    {
        Vector discretization_error = ZeroVector(rGaussPointsExact.size());

        double ele_error;

        for (IndexType tri_i = 0; tri_i < rGaussPointsExact.size(); ++tri_i)
        {
            ele_error = 0;
            for (IndexType point_i = 0; point_i < 3; ++point_i)
            {
                ele_error += sqrt(pow((rGaussPointsExact[tri_i](point_i, 0) - rGaussPointsApprox[tri_i](point_i, 0)),2) +
                                  pow((rGaussPointsExact[tri_i](point_i, 1) - rGaussPointsApprox[tri_i](point_i, 1)),2) +
                                  pow((rGaussPointsExact[tri_i](point_i, 2) - rGaussPointsApprox[tri_i](point_i, 2)),2));
            }
            discretization_error[tri_i] = ele_error;
        }
        return discretization_error;
    }



    /**
     * @brief This method initializes the necessary data structure for triangle
     * @param tr This is a struct used to pass data into and out of the triangulate() procedure
     */
    void CadTessellationModeler::InitTriangulationDataStructure(triangulateio& tr)
    {
        tr.pointlist                  = (REAL*) NULL;
        tr.pointattributelist         = (REAL*) NULL;
        tr.pointmarkerlist            = (int*) NULL;
        tr.numberofpoints             = 0;
        tr.numberofpointattributes    = 0;
        tr.trianglelist               = (int*) NULL;
        tr.triangleattributelist      = (REAL*) NULL;
        tr.trianglearealist           = (REAL*) NULL;
        tr.neighborlist               = (int*) NULL;
        tr.numberoftriangles          = 0;
        tr.numberofcorners            = 3;
        tr.numberoftriangleattributes = 0;
        tr.segmentlist                = (int*) NULL;
        tr.segmentmarkerlist          = (int*) NULL;
        tr.numberofsegments           = 0;
        tr.holelist                   = (REAL*) NULL;
        tr.numberofholes              = 0;
        tr.regionlist                 = (REAL*) NULL;
        tr.numberofregions            = 0;
        tr.edgelist                   = (int*) NULL;
        tr.edgemarkerlist             = (int*) NULL;
        tr.normlist                   = (REAL*) NULL;
        tr.numberofedges              = 0;
    };

    /**
     * @brief This method clears the data structure for triangle
     * @param tr This is a struct used to pass data into and out of the triangulate() procedure
     */
    void CadTessellationModeler::CleanTriangulationDataStructure(triangulateio& tr)
    {
        if(tr.pointlist != NULL)
        {
            free(tr.pointlist);
            tr.pointlist = nullptr;
        }
        if(tr.pointattributelist != NULL)
        {
            free(tr.pointattributelist);
            tr.pointattributelist = nullptr;
        }
        if(tr.pointmarkerlist != NULL)
        {
            free(tr.pointmarkerlist);
            tr.pointmarkerlist = nullptr;
        }
        if(tr.trianglelist != NULL)
        {
            free(tr.trianglelist);
            tr.trianglelist = nullptr;
        }
        if(tr.triangleattributelist != NULL)
        {
            free(tr.triangleattributelist);
            tr.triangleattributelist = nullptr;
        }
        if(tr.trianglearealist != NULL)
        {
            free(tr.trianglearealist);
            tr.trianglearealist = nullptr;
        }
        if(tr.neighborlist != NULL)
        {
            free(tr.neighborlist);
            tr.neighborlist = nullptr;
        }
        if(tr.segmentlist != NULL)
        {
            free(tr.segmentlist);
            tr.segmentlist = nullptr;
        }
        if(tr.segmentmarkerlist != NULL)
        {
            free(tr.segmentmarkerlist);
            tr.segmentmarkerlist = nullptr;
        }
        if(tr.holelist != NULL)
        {
            free(tr.holelist);
            tr.holelist = nullptr;
        }
        if(tr.regionlist != NULL)
        {
            free(tr.regionlist);
            tr.regionlist = nullptr;
        }
        if(tr.edgelist != NULL)
        {
            free(tr.edgelist);
            tr.edgelist = nullptr;
        }
        if(tr.edgemarkerlist != NULL)
        {
            free(tr.edgemarkerlist);
            tr.edgemarkerlist = nullptr;
        }
        if(tr.normlist != NULL)
        {
            free(tr.normlist);
            tr.normlist = nullptr;
        }
    }

    /**
     * @brief This method initializes the necessary data structure for triangle
     * @param tr This is a struct used to pass data into and out of the triangulate() procedure
     */
} // namespace Kratos
