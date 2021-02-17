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

        const auto& r_geometries = cad_model_part.Geometries();

        IndexType point_id = 0;
        std::vector<array_1d<double, 2>> boundary_loop_uv;

        for (auto it = r_geometries.begin(); it != r_geometries.end(); ++it) {
            if (it->GetGeometryType() == GeometryData::Kratos_Brep_Curve) {
                const auto& r_aux_geometry = *it;
                const auto r_brep_geom = dynamic_cast<const BrepCurveOnSurface<PointerVector<Node<3>>, PointerVector<Point>>&>(r_aux_geometry);
                const auto p_curve_on_surface_geom = r_brep_geom.pGetCurveOnSurface();

                auto tessellation = NurbsCurveTessellation<2, PointerVector<Node<3>>>::ComputeTessellation(
                    *p_curve_on_surface_geom,
                    p_curve_on_surface_geom->PolynomialDegree(0),
                    p_curve_on_surface_geom->DomainInterval(),
                    p_curve_on_surface_geom->KnotSpanIntervals(),
                    1e-2);

                boundary_loop_uv.resize(boundary_loop_uv.size() + tessellation.size() - 1, ZeroVector(2));
                for (int i = 0; i < tessellation.size() - 1; ++i)
                {
                    array_1d<double,3> result = ZeroVector(3);
                    array_1d<double,3> local_coordinate = ZeroVector(3);
                    local_coordinate[0] = tessellation.at(i).first;

                    for (IndexType j = 0; j < 2; ++j)
                    {
                        boundary_loop_uv[point_id][j] = p_curve_on_surface_geom->pGetCurve()->GlobalCoordinates(result, local_coordinate)[j];
                    }
                    point_id += 1;
                }
            }
        }

        // initiliazing the i/o containers
        struct triangulateio in_data;
        struct triangulateio out_data;
        struct triangulateio vor_out_data;

        InitTriangulationDataStructure(in_data);
        InitTriangulationDataStructure(out_data);
        InitTriangulationDataStructure(vor_out_data);

        // Initialize the pointlist (1d list) with the number of points and the coordinates
        // of the points (outer and inner polygons)
        in_data.numberofpoints = boundary_loop_uv.size();
        in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));
        in_data.pointmarkerlist = (int*) malloc(in_data.numberofpoints * sizeof(int));

        point_id = 0;
        IndexType point_marker_id = 0;
        IndexType point_marker = 0;

        for (IndexType node_i = 0; node_i < boundary_loop_uv.size(); ++node_i)
        {
            for (IndexType coords_i = 0; coords_i < 2; ++coords_i)
            {
                in_data.pointlist[point_id++] = boundary_loop_uv[node_i][coords_i];
            }
            in_data.pointmarkerlist[point_marker_id++] = point_marker;
        }

        // Initilize the segment list with the number of boundary edges and the start and end node id
        // For closed polygons the number of segments is equal to the number of points
        in_data.numberofsegments = boundary_loop_uv.size();
        in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
        in_data.segmentmarkerlist = (int*) malloc(in_data.numberofsegments * sizeof(int));
        IndexType node_id = 0;
        IndexType seg_marker = 0;
        IndexType start_node_id = 0;
        IndexType end_node_id = boundary_loop_uv.size();

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

        triangulate("Qqpza1", &in_data, &out_data, &vor_out_data);
        std::vector<Matrix> triangulation_uv;
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

        node_id = 0;
        IndexType vertex_id = 0;
        IndexType element_id = 0;

        for (int tri_i = 0; tri_i < triangulation_uv.size(); ++tri_i)
        {
            for (int point_i = 0; point_i < 3; ++point_i)
            {
                skin_model_part.CreateNewNode(node_id++,
                    triangulation_uv[tri_i](point_i,0),
                    triangulation_uv[tri_i](point_i,1),
                    0);
            }
        }

        Properties::Pointer p_properties(new Properties(0));

        // create elements in skin_model_part
        for (int element_i = 0; element_i < triangulation_uv.size(); ++element_i)
        {
            skin_model_part.CreateNewElement("Element3D3N",
                element_id++, {{vertex_id, vertex_id + 1, vertex_id + 2}}, p_properties);
            vertex_id += 3;
        }

    }


    /**
     * @brief This method returns the tessellation of the boundary curves
     */
    void CadTessellationModeler::ComputeBoundaryTessellation()
    {

    }

    /**
     * @brief This method returns the triangulation of a NURBS surface
     */
    void CadTessellationModeler::ComputeTriangulation()
    {

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
