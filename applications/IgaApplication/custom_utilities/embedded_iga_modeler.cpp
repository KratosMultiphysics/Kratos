//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
#include "embedded_iga_modeler.h"


namespace Kratos
{
    void EmbeddedIgaModeler::CreateElements2D(ModelPart& rSkinModelPart)
    {
        std::vector<array_1d<double,3>> polygon;
        EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector); 
        embedded_tessellation.CreateTessellationCurve(polygon);

        // Create Nodes in Skin Modelpart from the Points generated in the tessellation of the curve
        for (unsigned int point_i = 0; point_i < polygon.size(); ++point_i)
        {
            rSkinModelPart.CreateNewNode(point_i, polygon[point_i][0], polygon[point_i][1], polygon[point_i][2]);
        }

        // skin_model_part property Container (Container is empty, but for the skin no properties are needed)
        Properties::Pointer prop = rSkinModelPart.pGetProperties(0);

        unsigned int node_id = 0;
        // Create Elements in skin_model_part
        for (unsigned int element_i = 0; element_i < polygon.size() - 1; ++element_i)
        {
            rSkinModelPart.CreateNewElement("Element2D2N", element_i, {{node_id, node_id + 1}}, prop);
            node_id += 1;
        }
    }

    void EmbeddedIgaModeler::CreateElements3D(ModelPart& rSkinModelPart)
    {
        std::vector<array_1d<double,3>> triangulation_xyz; 
        MapTriangulationGeometricSpace(triangulation_xyz); 

        // for (unsigned int i = 0; i < triangulation_xyz.size(); ++i)     KRATOS_WATCH(triangulation_xyz[i])

        
        // Create Nodes in Skin Modelpart from the Points generated in the tessellation of the curve
        for (unsigned int point_i = 0; point_i < triangulation_xyz.size(); ++point_i)
        {
            rSkinModelPart.CreateNewNode(point_i, triangulation_xyz[point_i][0], triangulation_xyz[point_i][1], triangulation_xyz[point_i][2]);
        }

        // skin_model_part property Container (Container is empty, but for the skin no properties are needed)
        Properties::Pointer prop = rSkinModelPart.pGetProperties(0);

        unsigned int node_id = 0;
        const auto number_elements = triangulation_xyz.size() / 3; 

        KRATOS_WATCH(number_elements)
        // Create Elements in skin_model_part
        for (unsigned int element_i = 0; element_i < number_elements; ++element_i)
        {
            rSkinModelPart.CreateNewElement("Element3D3N", element_i, {{node_id, node_id + 1, node_id + 2}}, prop);
            node_id += 3;
        }

        KRATOS_WATCH(rSkinModelPart)
    }


    std::vector<Matrix> EmbeddedIgaModeler::TriangulateEmpire()
    {
        std::vector<array_1d<double,3>> polygon;
        EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector); 
        embedded_tessellation.CreateTessellationParameterCurve(polygon);
        
        std::vector<Matrix> triangles;
        EmbeddedIgaTriangulation embedded_triangulation; 
        embedded_triangulation.CreateTrianglesEmpire(polygon, triangles); 
        
        return triangles;
    }


    void EmbeddedIgaModeler::TriangulateCurveOnSurface(std::vector<array_1d<double,2>>& rTriangulation_uv)
    {
        /**
         * Triangulation of the patch in the parametric space
        */

        std::vector<array_1d<double,2>> polygon;
        EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector);
        embedded_tessellation.CreateTessellationCurveOnSurface(polygon);

        // initializing the i/o containers
        struct triangulateio in_data; 
        struct triangulateio out_data; 
        struct triangulateio vor_out_data;

        InitTriangulationDataStructure(in_data); 
        InitTriangulationDataStructure(out_data); 
        InitTriangulationDataStructure(vor_out_data); 

        // Initialize the pointlist (1d list) with the number of points and the position 
        in_data.numberofpoints = polygon.size(); 
        in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));

        unsigned int point_idx = 0;
        for (unsigned int i = 0; i < in_data.numberofpoints; ++i)
        {
            for (unsigned int j = 0; j < 2; ++j)    in_data.pointlist[point_idx++] = polygon[i][j];
        }
        
        // Initilize the segment list with the number of boundary edges and the start and end node id
        // For closed polygons the number of segments is equal to the number of points
        in_data.numberofsegments = in_data.numberofpoints; 
        in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
        
        unsigned int vertex_id = 1;
        for (unsigned int seg_idx = 1; seg_idx < in_data.numberofsegments * 2 - 1;)
        {
            for (unsigned int j = 0; j < 2; ++j)    in_data.segmentlist[seg_idx++] =  vertex_id;
            vertex_id += 1; 
        }
        in_data.segmentlist[0] = 0; 
        in_data.segmentlist[in_data.numberofsegments * 2 - 1] = 0; 

        char trigenOptsVerbose[] = "pz"; 
        char* trigenOpts = trigenOptsVerbose; 

        triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

        rTriangulation_uv.resize(out_data.numberoftriangles * 3); 

        KRATOS_WATCH(out_data.numberoftriangles)

        unsigned int tri_id = 0; 
        for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {
                rTriangulation_uv[tri_id + j][0] = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2];
                rTriangulation_uv[tri_id + j][1] = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2 + 1]; 
            }
            tri_id += 3;  
        }
    }

    void EmbeddedIgaModeler::MapTriangulationGeometricSpace(std::vector<array_1d<double,3>>& rTriangulation_xyz)
    {
        /**
         * Map Points generated in the triangulation of the patch from the 
         * parametric space into the geometric space
        */
        std::cout << "MapTriangulationGeometricSpace" << std::endl; 

        std::vector<array_1d<double,2>> triangulation_uv; 
        TriangulateCurveOnSurface(triangulation_uv);

        rTriangulation_xyz.resize(triangulation_uv.size());  

        for (unsigned int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
        {
            for (unsigned int face_i = 0; face_i < m_brep_model_vector[brep_i].GetFaceVector().size(); ++face_i)
            {
                auto geometry = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface(); 

                for (unsigned int point_i = 0; point_i < triangulation_uv.size(); ++point_i)
                {    
                    auto point_xyz = geometry->PointAt(triangulation_uv[point_i][0],triangulation_uv[point_i][1]);

                    for (unsigned int j = 0; j < 3; ++j)    rTriangulation_xyz[point_i][j] = point_xyz[j]; 
                }
            }
        }
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintTriangulationPoints_uv()
    {
        std::vector<array_1d<double,2>> triangulation_uv; 
        TriangulateCurveOnSurface(triangulation_uv); 

        std::vector<std::vector<double>> coords_uv(triangulation_uv.size(), std::vector<double>(2.0)); 
        
        for (unsigned int i = 0; i < triangulation_uv.size(); ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)    coords_uv[i][j] = triangulation_uv[i][j]; 
        }
        return coords_uv; 
    }


    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintTriangulationPoints_xyz()
    {
        std::vector<array_1d<double,3>> triangulation_xyz; 
        MapTriangulationGeometricSpace(triangulation_xyz);

        std::vector<std::vector<double>> coords_xyz(triangulation_xyz.size(), std::vector<double>(3,0)); 

        for (unsigned int i = 0; i < triangulation_xyz.size(); ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)    coords_xyz[i][j] = triangulation_xyz[i][j]; 
        }
        return coords_xyz; 
    }


    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintCurveTessellationPoints()
    {
        std::vector<array_1d<double,3>> polygon;
        EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector);
        embedded_tessellation.CreateTessellationCurve(polygon);
        std::vector<std::vector<double>> coords(polygon.size(), std::vector<double>(3,0));

        for (unsigned int i = 0; i < polygon.size(); ++i)
        {    
            for (unsigned int j = 0; j < 3; ++j)    coords[i][j] = polygon[i][j];
        }
        return coords; 
    }

    std::vector<std::vector<double>> EmbeddedIgaModeler::PrintCurveOnSurfaceTessellationPoints()
    {
        std::vector<array_1d<double,2>> polygon;
        EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector);
        embedded_tessellation.CreateTessellationCurveOnSurface(polygon);
        std::vector<std::vector<double>> coords(polygon.size(), std::vector<double>(3,0));

        for (unsigned int i = 0; i < polygon.size(); ++i)
        {    
            for (unsigned int j = 0; j < 2; ++j)    coords[i][j] = polygon[i][j];
        }
        return coords; 
    }



    // std::vector<std::vector<double>> EmbeddedIgaModeler::PrintParameterCurveTessellationPoints()
    // {
    //     std::vector<array_1d<double,3>> polygon;
    //     EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector);
    //     embedded_tessellation.CreateTessellationParameterCurve(polygon);
        
    //     std::vector<std::vector<double>> coords(polygon.size(), std::vector<double>(2,0));

    //     for (unsigned int i = 0; i < polygon.size(); ++i)
    //     {
    //         for (unsigned int j = 0; j < 2; ++j)    coords[i][j] = polygon[i][j];
    //     }
    //     return coords;
    // }

    // std::vector<std::vector<double>> EmbeddedIgaModeler::PrintGaussPoints()
    // {
    //     std::vector<Matrix> triangles = TriangulateEmpire(); 
    //     EmbeddedIgaErrorEstimation error_estimator(triangles); 
    //     std::vector<array_1d<double, 2> > gp_pos; 
    //     error_estimator.InsertGaussPoints(gp_pos);

    //     std::vector<std::vector<double>> gp_coords(gp_pos.size(), std::vector<double>(2,0)); 
        
    //     for (unsigned int i = 0; i < gp_pos.size(); ++i)
    //     {
    //         for (unsigned int j = 0; j < 2; ++j)    gp_coords[i][j] = gp_pos[i][j];
    //     }
    //     return gp_coords; 
    // }

    // std::vector<std::vector<double>> EmbeddedIgaModeler::PrintMappedGaussPoints()
    // {
    //     std::vector<Matrix> triangles = TriangulateEmpire(); 
    //     EmbeddedIgaErrorEstimation error_estimator(triangles); 
    //     std::vector<array_1d<double, 2> > gp_pos; 

    //     error_estimator.InsertGaussPoints(gp_pos); 
        
    //     auto geometry = m_brep_model_vector[0].GetFaceVector()[0].GetSurface(); 

    //     std::vector<std::vector<double>> mapped_gp_coords (gp_pos.size(), std::vector<double>(3,0)); 

    //     for (unsigned int i = 0; i < gp_pos.size(); ++i)
    //     {    
    //         auto point = geometry->PointAt(gp_pos[i][0],gp_pos[i][1]); 
    //         for (unsigned int j = 0; j < 3; ++j)    mapped_gp_coords[i][j] = point[j]; 
    //     }
    //     return mapped_gp_coords; 
    // }

    

    // std::vector<std::vector<double>> EmbeddedIgaModeler::TestTriangle()
    // {
    //     std::vector<array_1d<double,2>> polygon;

    //     polygon.resize(5); 

    //     polygon[0][0] = 0;
    //     polygon[0][1] = 0; 
    //     polygon[1][0] = 10;
    //     polygon[1][1] = 0; 
    //     polygon[2][0] = 10;
    //     polygon[2][1] = 5; 
    //     polygon[3][0] = 0;
    //     polygon[3][1] = 5;
    //     polygon[4][0] = 7; 
    //     polygon[4][1] = 3; 




    //     // initializing the i/o containers
    //     struct triangulateio in_data; 
    //     struct triangulateio out_data; 
    //     struct triangulateio vor_out_data;

    //     vor_out_data.numberofpoints = 1; 
    //     vor_out_data.pointlist = (REAL*) malloc(vor_out_data.numberofpoints * 2 * sizeof(REAL));
    //     vor_out_data.pointlist[0] = 8; 
    //     vor_out_data.pointlist[1] = 3; 

    //     InitTriangulationDataStructure(in_data); 
    //     InitTriangulationDataStructure(out_data); 
    //     InitTriangulationDataStructure(vor_out_data); 

    //     // Initialize the pointlist (1d list) with the number of points and the position 
    //     in_data.numberofpoints = polygon.size(); 
    //     in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));

    //     unsigned int point_idx = 0;
    //     for (unsigned int i = 0; i < in_data.numberofpoints; ++i)
    //     {
    //         for (unsigned int j = 0; j < 2; ++j)    in_data.pointlist[point_idx++] = polygon[i][j];
    //     }
        
    //     // Initilize the segment list with the number of boundary edges and the start and end node id
    //     // For closed polygons the number of segments is equal to the number of points
    //     in_data.numberofsegments = in_data.numberofpoints; 
    //     in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
        
    //     unsigned int vertex_id = 1;
    //     for (unsigned int seg_idx = 1; seg_idx < in_data.numberofsegments * 2 - 1;)
    //     {
    //         for (unsigned int j = 0; j < 2; ++j)    in_data.segmentlist[seg_idx++] =  vertex_id;
    //         vertex_id += 1; 
    //     }
    //     in_data.segmentlist[0] = 0; 
    //     in_data.segmentlist[in_data.numberofsegments * 2 - 1] = 0; 

    //     char trigenOptsVerbose[] = "Vvpz"; 
    //     char* trigenOpts = trigenOptsVerbose; 

    //     triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

    //     std::vector<std::vector<double>> tri_coords (out_data.numberoftriangles * 3 , std::vector<double>(2,0)); 
    //     unsigned int id = 0; 
    //     for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
    //     {
    //         for (unsigned int j = 0; j < 3; ++j)
    //         {
    //             tri_coords[id + j][0] = out_data.pointlist[out_data.trianglelist[id + j] * 2];
    //             tri_coords[id + j][1] = out_data.pointlist[out_data.trianglelist[id + j] * 2 + 1]; 
    //         }
    //         id += 3; 
    //     }
    //     return tri_coords; 
    // }

    // std::vector<std::vector<double>> EmbeddedIgaModeler::Triangulate()
    // {
    //     std::vector<array_1d<double,3>> polygon;
    //     EmbeddedIgaTessellation embedded_tessellation(m_brep_model_vector);
    //     embedded_tessellation.CreateTessellationCurve(polygon);

    //     // initializing the i/o containers
    //     struct triangulateio in_data; 
    //     struct triangulateio out_data; 
    //     struct triangulateio vor_out_data;

    //     InitTriangulationDataStructure(in_data); 
    //     InitTriangulationDataStructure(out_data); 
    //     InitTriangulationDataStructure(vor_out_data); 

    //     // Initialize the pointlist (1d list) with the number of points and the position 
    //     in_data.numberofpoints = polygon.size(); 
    //     in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));

    //     unsigned int point_idx = 0;
    //     for (unsigned int i = 0; i < in_data.numberofpoints; ++i)
    //     {
    //         for (unsigned int j = 0; j < 2; ++j)    in_data.pointlist[point_idx++] = polygon[i][j];
    //     }
        
    //     // Initilize the segment list with the number of boundary edges and the start and end node id
    //     // For closed polygons the number of segments is equal to the number of points
    //     in_data.numberofsegments = in_data.numberofpoints; 
    //     in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
        
    //     unsigned int vertex_id = 1;
    //     for (unsigned int seg_idx = 1; seg_idx < in_data.numberofsegments * 2 - 1;)
    //     {
    //         for (unsigned int j = 0; j < 2; ++j)    in_data.segmentlist[seg_idx++] =  vertex_id;
    //         vertex_id += 1; 
    //     }
    //     in_data.segmentlist[0] = 0; 
    //     in_data.segmentlist[in_data.numberofsegments * 2 - 1] = 0; 

    //     char trigenOptsVerbose[] = "pz"; 
    //     char* trigenOpts = trigenOptsVerbose; 

    //     triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

    //     std::vector<std::vector<double>> tri_coords (out_data.numberoftriangles * 3 , std::vector<double>(2,0)); 
    //     unsigned int id = 0; 
    //     for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
    //     {
    //         for (unsigned int j = 0; j < 3; ++j)
    //         {
    //             tri_coords[id + j][0] = out_data.pointlist[out_data.trianglelist[id + j] * 2];
    //             tri_coords[id + j][1] = out_data.pointlist[out_data.trianglelist[id + j] * 2 + 1]; 
    //         }
    //         id += 3; 
    //     }
    //     return tri_coords; 
    // }

    

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
        : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
        m_model_part(rModelPart)
    {
    }
} // namespace Kratos.