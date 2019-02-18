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
#include "embedded_iga_tessellation.h"

namespace Kratos
{
std::vector<std::vector<double>> EmbeddedIgaTessellation::CreateTessellation(
    const BrepFace& rFaceGeometry)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// TESSELLATION
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<array_1d<double, 2>> outer_polygon;
    std::vector<std::vector<array_1d<double, 2>>> inner_polygon;
     



    const auto surface_geometry = rFaceGeometry.GetSurface();
    for (unsigned int b_loop_i = 0; b_loop_i < rFaceGeometry.GetBoundaryLoop().size(); ++b_loop_i)
    {
        unsigned int point_index = 0;
        auto boundary_loop = rFaceGeometry.GetBoundaryLoop()[b_loop_i];

        std::vector<array_1d<double, 2>> inner_loop; //


        for (unsigned int trim_i = 0; trim_i < boundary_loop.GetTrimmingCurves().size(); ++trim_i)
        {
            const auto trimming_curve = boundary_loop.GetTrimmingCurves()[trim_i];
            int trim_index = trimming_curve.GetTrimIndex();

            const auto curve_geometry = rFaceGeometry.GetTrimCurve(trim_index);

            const auto curve_on_surface = CurveOnSurface<3>(
                curve_geometry->CurveGeometry(),
                surface_geometry,
                curve_geometry->Domain());

            const auto tessellation = Kratos::make_shared<ANurbs::CurveTessellation<array_1d<double, 3>>>();

            tessellation->Compute(curve_on_surface, 1e-2);

            if (boundary_loop.IsOuterLoop())
            {
                // polygon vector needs to be resized to the current number of points + the new points
                // generated in the tessellation. However, the first point of one trimming curve
                // is the last point of the previous trimming curve. To account for these points
                // we subtract - 1.
                outer_polygon.resize(outer_polygon.size() + tessellation->NbPoints() - 1);
                
                for (unsigned int i = 0; i < tessellation->NbPoints() - 1; ++i)
                {
                    for (unsigned int j = 0; j < 2; ++j)
                    {
                        outer_polygon[point_index][j] = curve_geometry->CurveGeometry()
                                                           ->PointAt(tessellation->Parameter(i))[j];
                    }
                    point_index += 1;
                }
            }
            else 
            {
            // if boundary_loop.IsOuterLoop() is false the loop is an inner loop and needs to 
            // be considered in the triangulation as a hole
                inner_loop.resize(inner_loop.size() + tessellation->NbPoints() - 1); 
                for (unsigned int i = 0; i < tessellation->NbPoints() - 1; ++i)
                {
                    for (unsigned int j = 0; j < 2; ++j)
                    {
                        inner_loop[point_index][j] = curve_geometry->CurveGeometry()
                                                           ->PointAt(tessellation->Parameter(i))[j];
                    }
                    point_index += 1;
                }
            }
        }
        
        outer_polygon.resize(4); 

        outer_polygon[0][0] = 0;
        outer_polygon[0][1] = 0;
        outer_polygon[1][0] = 10;
        outer_polygon[1][1] = 0;
        outer_polygon[2][0] = 10;
        outer_polygon[2][1] = 5;
        outer_polygon[3][0] = 0;
        outer_polygon[3][1] = 5;

        
        inner_loop.resize(3);
        inner_loop[0][0] = 2;
        inner_loop[0][1] = 2; 
        inner_loop[1][0] = 4;
        inner_loop[1][1] = 2;
        inner_loop[2][0] = 3;
        inner_loop[2][1] = 3;  





        inner_polygon.push_back(inner_loop); 
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// TRIANGULATION
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /**
     * This function generates a triangulation of the patch in the parametric space
    */


    // initializing the i/o containers
    struct triangulateio in_data; 
    struct triangulateio out_data; 
    struct triangulateio vor_out_data;
    InitTriangulationDataStructure(in_data); 
    InitTriangulationDataStructure(out_data); 
    InitTriangulationDataStructure(vor_out_data); 

    
    // Initialize the pointlist (1d list) with the number of points and the coordinates
    // of the points (outer and inner polygons) 
    
    unsigned int number_points = 0; 
    for (unsigned int i = 0; i < inner_polygon.size(); ++i)
    {
        number_points += inner_polygon[i].size(); 
    }
    number_points += outer_polygon.size(); 

    in_data.numberofpoints = number_points; 
    in_data.pointlist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));
    in_data.pointmarkerlist = (int*) malloc(in_data.numberofpoints * sizeof(int));

    unsigned int point_index = 0;
    unsigned int point_marker = 0; 
    for (unsigned int node_i = 0; node_i < outer_polygon.size(); ++node_i)
    {
        for (unsigned int coords_i = 0; coords_i < 2; ++coords_i)    
        {
            in_data.pointlist[point_index] = outer_polygon[node_i][coords_i];
            in_data.pointmarkerlist[point_index] = point_marker; 
            point_index += 1; 
        }
    }

    for (unsigned int poly_i; poly_i < inner_polygon.size(); ++poly_i)
    {
        point_marker += 1; 
        for (unsigned int node_i = 0; node_i < inner_polygon[poly_i].size(); ++node_i)
        {
            for (unsigned int coords_i = 0; coords_i < 2; ++coords_i)    
            {
                in_data.pointlist[point_index] = inner_polygon[poly_i][node_i][coords_i];
                in_data.pointmarkerlist[point_index] = point_marker; 
                point_index += 1; 
            }
        }
    }

    // Initilize the segment list with the number of boundary edges and the start and end node id
    // For closed polygons the number of segments is equal to the number of points
    in_data.numberofsegments = number_points; 
    in_data.segmentlist = (int*) malloc(in_data.numberofsegments * 2 * sizeof(int));
    in_data.segmentmarkerlist = (int*) malloc(in_data.numberofsegments * sizeof(int));
    
    unsigned int node_id = 0; 
    unsigned int seg_marker = 0; 
    int seg_i = 0; 
    for (unsigned int i = 0; i < outer_polygon.size() * 2; ++i)
    {
        for (unsigned int j = 0; j < 2; ++j)
        {
            in_data.segmentlist[seg_i] = node_id; 
            in_data.segmentmarkerlist[seg_i] = seg_marker;
            
            if ((seg_i == 0)) 
            {
                break; 
            }
            if (seg_i == outer_polygon.size() * 2 - 1)
            {
                in_data.segmentlist[outer_polygon.size() - 1] = 0;
                break; 
            }
            KRATOS_WATCH(seg_i)
            seg_i += 1; 
        }
        node_id += 1; 
    }

    for(unsigned int i = 0; i < in_data.numberofsegments * 2; ++i)
    {
        KRATOS_WATCH(in_data.segmentlist[i])
    }


    for(unsigned int i = 0; i < in_data.numberofsegments * 2; ++i)
    {
        KRATOS_WATCH(in_data.segmentmarkerlist[i])
    }


    // unsigned int node_id = 0;
    // unsigned int seg_marker = 0;  
    
    // in_data.segmentlist[0] = node_id; 
    // in_data.segmentmarkerlist[0] = seg_marker;

    // for (unsigned int seg_i = 0; seg_i < outer_polygon.size() * 2;)
    // {
    //     node_id += 1; 
    //     for (unsigned int j = 0; j < 2; ++j)
    //     {
    //         if (seg_i < outer_polygon.size() * 2 - 1)
    //         {
    //             in_data.segmentlist[seg_i] = node_id; 
    //         }
    //         else
    //         {
    //             in_data.segmentlist[seg_i] = in_data.segmentlist[0]; 
    //         }
            
    //         in_data.segmentmarkerlist[seg_i] = seg_marker; 
    //         seg_i += 1; 
    //     }
    // }


    // for (unsigned int poly_i = 0; poly_i < inner_polygon.size(); ++poly_i)
    // {        
    //     node_id += 1;  
    //     seg_marker += 1;
    //     unsigned int start_id = node_id;

    //     in_data.segmentlist[seg_i] = node_id; 
    //     in_data.segmentmarkerlist[0] = seg_marker; 

    //     for (unsigned int inner_seg_i = seg_i; seg_i < )
    // }







    // unsigned int node_id = 1; 
    // unsigned int seg_i;
    // unsigned int seg_marker = 0;  
    // in_data.segmentlist[0] = 0;
    // in_data.segmentmarkerlist[0] = seg_marker;
    // for (seg_i = 1; seg_i < outer_polygon.size() * 2 - 1;)
    // {
    //     for (unsigned int j = 0; j < 2; ++j)
    //     {
    //         in_data.segmentlist[seg_i] = node_id; 
    //         in_data.segmentmarkerlist[seg_i] = seg_marker; 
    //         seg_i += 1;
    //     }
    //     node_id += 1; 
    // }
    // in_data.segmentlist[seg_i] = 0;
    // in_data.segmentmarkerlist[seg_i++] = seg_marker;
    

    // for (unsigned int poly_i = 0; poly_i < inner_polygon.size(); ++poly_i)
    // {
    //     seg_marker += 1; 
    //     for (inner_seg_i = seg_i; inner_seg_i < )
    // }
    






    in_data.numberofholes = inner_polygon.size(); 
    in_data.holelist = (REAL*) malloc(in_data.numberofpoints * 2 * sizeof(REAL));




    char trigenOptsVerbose[] = "Dpza1"; 
    char* trigenOpts = trigenOptsVerbose; 

    triangulate(trigenOpts, &in_data, &out_data, &vor_out_data);

    std::vector<Matrix> triangulation(
        out_data.numberoftriangles, ZeroMatrix(3,2)); 

    unsigned int tri_id = 0; 
    for (unsigned int i = 0; i < out_data.numberoftriangles; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {   
            triangulation[i](j,0) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2];
            triangulation[i](j,1) = out_data.pointlist[out_data.trianglelist[tri_id + j] * 2 + 1]; 
        }
        tri_id += 3;  
    }
    
    std::vector<std::vector<double>> coords_uv(triangulation.size() * 3, std::vector<double>(2,0)); 
        
    unsigned int point_i = 0; 
    for (unsigned int i = 0; i < triangulation.size(); ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)    
        {
            for (unsigned int k = 0; k < 2; ++k)    coords_uv[point_i][k] = triangulation[i](j,k); 
            point_i += 1; 
        }
    }
    return coords_uv; 

}

EmbeddedIgaTessellation::EmbeddedIgaTessellation()

{
}
} // namespace Kratos.