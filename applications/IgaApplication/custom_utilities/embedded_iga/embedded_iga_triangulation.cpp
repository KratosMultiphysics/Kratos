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
#include "embedded_iga_triangulation.h"

namespace Kratos
{
void EmbeddedIgaTriangulation::CreateTriangulation(
    const double mTriangulationError,
    const double mInitialTriangleArea,
    const int mMaxTriangulationIterations,
    const int mEchoLevel,
    const BrepFace& rFaceGeometry,
    const std::vector<std::vector<array_1d<double,2>>>& rOuterPolygon,
    const std::vector<std::vector<array_1d<double,2>>>& rInnerPolygon,
    std::vector<Matrix>& rTriangulation_xyz)
{
    /**
     * This function generates a triangulation of the patch in the parametric space
    */

    // initializing the i/o containers
    struct triangulateio in_data; 
    struct triangulateio out_data; 
    struct triangulateio vor_out_data;

    InitTriangulationDataStructure(in_data); 

    // Initialize the pointlist (1d list) with the number of points and the coordinates
    // of the points (outer and inner polygons) 
    unsigned int number_points = 0; 
    for (unsigned int i = 0; i < rOuterPolygon.size(); ++i)
    {
        number_points += rOuterPolygon[i].size(); 
    }
    for (unsigned int i = 0; i < rInnerPolygon.size(); ++i)
    {
        number_points += rInnerPolygon[i].size(); 
    }
    
    in_data.numberofpoints = number_points; 
    in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));
    in_data.pointmarkerlist = (int*) malloc(in_data.numberofpoints * sizeof(int));

    unsigned int point_index = 0;
    unsigned int point_marker_index = 0;
    unsigned int point_marker = 0; 
    for (unsigned int poly_i = 0; poly_i < rOuterPolygon.size(); ++poly_i)
    {
        for (unsigned int node_i = 0; node_i < rOuterPolygon[poly_i].size(); ++node_i)    
        {
            for (unsigned int coords_i = 0; coords_i < 2; ++coords_i)
            {
                in_data.pointlist[point_index++] = rOuterPolygon[poly_i][node_i][coords_i];
            }
            in_data.pointmarkerlist[point_marker_index++] = point_marker; 
        }
        point_marker++;
    }

    for (unsigned int poly_i = 0; poly_i < rInnerPolygon.size(); ++poly_i)
    {
        for (unsigned int node_i = 0; node_i < rInnerPolygon[poly_i].size(); ++node_i)
        {
            for (unsigned int coords_i = 0; coords_i < 2; ++coords_i)    
            {
                in_data.pointlist[point_index++] = rInnerPolygon[poly_i][node_i][coords_i]; 
            }
            in_data.pointmarkerlist[point_marker_index++] = point_marker; 
        }
        point_marker++;
    }

    // Initilize the segment list with the number of boundary edges and the start and end node id
    // For closed polygons the number of segments is equal to the number of points
    in_data.numberofsegments = number_points; 
    in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
    in_data.segmentmarkerlist = (int*) malloc(in_data.numberofsegments * sizeof(int));
    
    unsigned int node_id = 0; 
    unsigned int seg_marker = 0; 
    unsigned int start_node_id = 0;
    unsigned int end_node_id = 0;
     
    for (unsigned int poly_i = 0; poly_i < rOuterPolygon.size(); ++poly_i)
    {
        end_node_id += rOuterPolygon[poly_i].size(); 

        for (unsigned int seg_i = start_node_id * 2 ; seg_i < end_node_id * 2; ++seg_i)
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
        seg_marker++;
        start_node_id = end_node_id;
    }

    for (unsigned int poly_i = 0; poly_i < rInnerPolygon.size(); ++poly_i)
    {
        end_node_id += rInnerPolygon[poly_i].size(); 

        for (unsigned int seg_i = start_node_id * 2 ; seg_i < end_node_id * 2; ++seg_i)
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
        seg_marker++;
        start_node_id = end_node_id;
    }

    in_data.numberofholes = rInnerPolygon.size(); 
    in_data.holelist = (REAL*) malloc(in_data.numberofholes * 2 * sizeof(REAL));

    for (unsigned int poly_i = 0; poly_i < rInnerPolygon.size(); ++poly_i)
    {
        
        
        mapbox::geometry::polygon<double> polygon; 
        mapbox::geometry::linear_ring<double> linear_ring; 
        for (unsigned int node_i = 0; node_i < rInnerPolygon[poly_i].size(); ++node_i)
        {
            linear_ring.push_back(mapbox::geometry::point<double> (rInnerPolygon[poly_i][node_i][0], rInnerPolygon[poly_i][node_i][1]));
        }
        polygon.push_back(linear_ring); 
     
        mapbox::geometry::point<double> inner_coord = mapbox::polylabel(polygon, 1e-6);
        in_data.holelist[poly_i * 2] = inner_coord.x;
        in_data.holelist[poly_i * 2 + 1] = inner_coord.y;
    }

    double max_area = mInitialTriangleArea;
    
    std::vector<Matrix> triangulation_uv; 
    for (int it = 1; it < mMaxTriangulationIterations; it++)
    {
        InitTriangulationDataStructure(out_data); 
        InitTriangulationDataStructure(vor_out_data); 
        
        
        std::ostringstream stream_obj;
        stream_obj << std::fixed;
        stream_obj << std::setprecision(16);
        stream_obj << max_area;
        std::string area_str = stream_obj.str();  
        int n = area_str.length(); 
        char char_array[n + 1];      
        strcpy(char_array, area_str.c_str());
        char trigenOpts[25] = "Qqpza";
        strcat(trigenOpts, char_array); 
        
        triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

        triangulation_uv.resize(out_data.numberoftriangles, ZeroMatrix(3,2));

        unsigned int tri_id = 0; 
        for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
            {   
                triangulation_uv[i](j,0) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2];
                triangulation_uv[i](j,1) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2 + 1]; 
            }
            tri_id += 3;  
        }

        std::vector<Matrix> gauss_points_exact_xyz; 
        EmbeddedIgaErrorEstimation::InsertGaussPointsExactSurface(
            rFaceGeometry, triangulation_uv, gauss_points_exact_xyz);

        std::vector<Matrix> gauss_points_approx_xyz; 
        EmbeddedIgaErrorEstimation::InsertGaussPointsApproxSurface(
            rFaceGeometry, triangulation_uv, gauss_points_approx_xyz); 
            
        Vector error; 
        EmbeddedIgaErrorEstimation::GetError(
            gauss_points_exact_xyz, gauss_points_approx_xyz, error);
        
        auto max_error = *std::max_element(std::begin(error), std::end(error));
        
        auto tolerance = false; 
        for (unsigned int i = 0; i < error.size(); ++i)
        {
            if (error[i] > mTriangulationError)
            {
                tolerance = true; 
                break; 
            }
        }
        KRATOS_INFO_IF("EMBEDDED_IGA", mEchoLevel >= 1) << "Iteration " << it << std::endl;
        KRATOS_INFO_IF("EMBEDDED_IGA", mEchoLevel >= 1) << "Area: " << max_area << " - max_error: " << max_error << std::endl;
        

        // Triangle copies the pointer for the holelist from the in_data to the out_data
        // In order to avoid the freeing of memory twice, which leads to an error the points will
        // be deleted from the out_data and in_data cleans the memory fully. 
        out_data.holelist = nullptr;
        
        CleanTriangulationDataStructure(out_data); 
        CleanTriangulationDataStructure(vor_out_data); 
        
        if (tolerance == false )    break; 
        max_area /= 2; 
    }
    CleanTriangulationDataStructure(in_data); 
    
    EmbeddedIgaMapper::MapCartesianSpace(
        rFaceGeometry, triangulation_uv, rTriangulation_xyz); 
}

EmbeddedIgaTriangulation::EmbeddedIgaTriangulation()
{}

} // namespace Kratos.