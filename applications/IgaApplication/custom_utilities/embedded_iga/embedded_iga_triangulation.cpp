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
    const std::vector<std::vector<array_1d<double,2>>>& rOuterPolygon,
    const std::vector<std::vector<array_1d<double,2>>>& rInnerPolygon)
{
    /**
     * This function generates a triangulation of the patch in the parametric space
    */

    // // initializing the i/o containers
    // struct triangulateio in_data; 
    // struct triangulateio out_data; 
    // struct triangulateio vor_out_data;

    // InitTriangulationDataStructure(in_data); 
    // InitTriangulationDataStructure(out_data); 
    // InitTriangulationDataStructure(vor_out_data); 

    // // Initialize the pointlist (1d list) with the number of points and the coordinates
    // // of the points (outer and inner polygons) 
    // unsigned int number_points = 0; 
    // for (unsigned int i = 0; i < rOuterPolygon.size(); ++i)
    // {
    //     number_points += rOuterPolygon[i].size(); 
    // }
    // for (unsigned int i = 0; i < rInnerPolygon.size(); ++i)
    // {
    //     number_points += rInnerPolygon[i].size(); 
    // }
    
    // in_data.numberofpoints = number_points; 
    // in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));
    // in_data.pointmarkerlist = (int*) malloc(in_data.numberofpoints * sizeof(int));

    // KRATOS_WATCH(in_data.numberofpoints)


    // unsigned int point_index = 0;
    // unsigned int point_marker_index = 0;
    // unsigned int point_marker = 0; 
    // for (unsigned int node_i = 0; node_i < outer_polygon.size(); ++node_i)
    // {
    //     for (unsigned int coords_i = 0; coords_i < 2; ++coords_i)    
    //     {
    //         in_data.pointlist[point_index] = outer_polygon[node_i][coords_i];
    //         point_index += 1; 
    //     }
    //     in_data.pointmarkerlist[point_marker_index] = point_marker; 
    //     point_marker_index++;
    // }
    // for (unsigned int poly_i = 0; poly_i < inner_polygon.size(); ++poly_i)
    // {
    //     point_marker += 1; 
    //     for (unsigned int node_i = 0; node_i < inner_polygon[poly_i].size(); ++node_i)
    //     {
    //         for (unsigned int coords_i = 0; coords_i < 2; ++coords_i)    
    //         {
    //             in_data.pointlist[point_index] = inner_polygon[poly_i][node_i][coords_i];
    //             point_index += 1; 
    //         }
    //         in_data.pointmarkerlist[point_marker_index] = point_marker; 
    //         point_marker_index++;
    //     }
    // }

    // // Initilize the segment list with the number of boundary edges and the start and end node id
    // // For closed polygons the number of segments is equal to the number of points
    // in_data.numberofsegments = number_points; 
    // in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
    // in_data.segmentmarkerlist = (int*) malloc(in_data.numberofsegments * sizeof(int));
    
    // unsigned int seg_index = 0; 
    // unsigned int node_id = 0; 
    // unsigned int seg_marker = 0; 
    // unsigned int end_node_id = outer_polygon.size(); 
    
    // for (unsigned int i = 0; i < outer_polygon.size() * 2; ++i)
    // {
    //     in_data.segmentlist[i] = node_id; 
        
    //     if (node_id == end_node_id)    
    //     {
    //         in_data.segmentlist[i] = 0; 
    //     }

    //     if (i % 2 == 0)   
    //     {
    //         in_data.segmentmarkerlist[i/2] = seg_marker; 
    //         node_id += 1; 
    //     }
    // }
    
    // number_points = 0;
    // unsigned int start_node_id = 0;
    // end_node_id = 0;
    
    // for (unsigned int i = 0; i < inner_polygon.size(); ++i)
    // {
    //     seg_marker += 1; // increase the marker for each polygon by one
    //     if (i != 0)   
    //     {
    //         number_points = 0; 
    //         for (unsigned int j = 0; j < i; ++j)
    //         {
    //             number_points += inner_polygon[j].size();  // sum up points of the inner hole before to find the start_node_id
    //         }
    //     }

    //     start_node_id = outer_polygon.size() + number_points;  
    //     end_node_id = start_node_id + inner_polygon[i].size(); 
    //     node_id = start_node_id; 

    //     for (unsigned int i = start_node_id * 2; i < end_node_id * 2; ++i)
    //     {
    //         in_data.segmentlist[i] = node_id;
            
    //         if (node_id == end_node_id)    
    //         {   
    //             in_data.segmentlist[i] = start_node_id; 
    //         }
    //         if (i % 2 == 0)   
    //         {
    //             in_data.segmentmarkerlist[i/2] = seg_marker; 
    //             node_id += 1; 
    //         }
    //     }
    // }

    // in_data.numberofholes = inner_polygon.size(); 
    // in_data.holelist = (REAL*) malloc(in_data.numberofholes * 2 * sizeof(REAL));

    // in_data.holelist[0] = 3; 
    // in_data.holelist[1] = 2.5;
    // in_data.holelist[2] = 6; 
    // in_data.holelist[3] = 3; 

    // KRATOS_WATCH(in_data.numberofpoints)
    // KRATOS_WATCH(in_data.numberofsegments)
    // KRATOS_WATCH(in_data.numberofholes)

    // for (int i = 0; i < in_data.numberofpoints * 2; ++i)    KRATOS_WATCH(in_data.pointlist[i])
    // for (int i = 0; i < in_data.numberofpoints; ++i)    KRATOS_WATCH(in_data.pointmarkerlist[i])
    // for (int i = 0; i < in_data.numberofsegments * 2; ++i)    KRATOS_WATCH(in_data.segmentlist[i])
    // for (int i = 0; i < in_data.numberofsegments; ++i)    KRATOS_WATCH(in_data.segmentmarkerlist[i])
    // for (int i = 0; i < in_data.numberofholes * 2; ++i)    KRATOS_WATCH(in_data.holelist[i])
    





    
    // char trigenOptsVerbose[] = "Dpza1"; 
    // char* trigenOpts = trigenOptsVerbose; 

    // triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

    // triangulation.resize(out_data.numberoftriangles, ZeroMatrix(3,2)); 

    // unsigned int tri_id = 0; 
    // for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
    // {
    //     for (unsigned int j = 0; j < 3; ++j)
    //     {   
    //         triangulation[i](j,0) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2];
    //         triangulation[i](j,1) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2 + 1]; 
    //     }
    //     tri_id += 3;  
    // }
    
    // KRATOS_WATCH(triangulation)
    // // std::vector<std::vector<double>> coords_uv;

  

    
    // // CleanTriangulationDataStructure(in_data); 
    // // CleanTriangulationDataStructure(out_data); 
    // // CleanTriangulationDataStructure(vor_out_data); 
    
    // // return coords_uv;
    
}

EmbeddedIgaTriangulation::EmbeddedIgaTriangulation()
{}

} // namespace Kratos.