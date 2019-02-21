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
void EmbeddedIgaModeler::Test()
{
    const auto face = m_brep_model_vector[0].GetFaceVector()[0];
    std::vector<std::vector<array_1d<double, 2>>> outer_polygon_uv;
    std::vector<std::vector<array_1d<double, 2>>> inner_polygon_uv;
    
    EmbeddedIgaTessellation::CreateTessellation(
        face, outer_polygon_uv, inner_polygon_uv); 

    std::vector<Matrix> triangulation_uv;

    EmbeddedIgaTriangulation embedded_triangulation; 
    embedded_triangulation.CreateTriangulation(
        outer_polygon_uv, inner_polygon_uv, triangulation_uv);

    std::vector<Matrix> gauss_points_uv; 
    EmbeddedIgaErrorEstimation::InsertGaussPointsExactSurface(
        triangulation_uv, gauss_points_uv); 

    // project the Gauss-Points onto the exact surface using the MapCartesianSpace member function
    std::vector<Matrix> gauss_points_xyz; 
    MapCartesianSpace(face, gauss_points_uv, gauss_points_xyz);
    
    KRATOS_WATCH(gauss_points_uv)
    KRATOS_WATCH(gauss_points_xyz)
    
}



void EmbeddedIgaModeler::MapCartesianSpace(
    const BrepFace& rFaceGeometry,
    const std::vector<Matrix>& rPoints_uv,
    std::vector<Matrix>& rPoints_xyz) 
{
    /**
     * This function maps points from the parametric space into the cartesian space.
     * The number of points per row (element of the std::vector) is dependent on the
     * number of input points.
    */
    const auto number_triangles = rPoints_uv.size(); 
    const auto surface_geometry = rFaceGeometry.GetSurface();
    
    rPoints_xyz.resize(
        rPoints_uv.size(), ZeroMatrix(rPoints_uv[0].size1(),3));


    for (unsigned int tri_i = 0; tri_i < number_triangles; ++tri_i)
    {    
        for (unsigned int point_i = 0; point_i < 3; ++point_i)
        {
            auto point_xyz = surface_geometry->PointAt(
                    rPoints_uv[tri_i](point_i,0), rPoints_uv[tri_i](point_i,1));
            
            for (unsigned int j = 0; j < 3; ++j)    
            {
                rPoints_xyz[tri_i](point_i, j) = point_xyz[j]; 
            }
        }
    }      
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<double>> EmbeddedIgaModeler::PrintParametricTessellation()
{
    const auto face = m_brep_model_vector[0].GetFaceVector()[0];
    
    std::vector<std::vector<array_1d<double, 2>>> outer_polygon;
    std::vector<std::vector<array_1d<double, 2>>> inner_polygon;
    
    EmbeddedIgaTessellation::CreateTessellation(
        face, outer_polygon, inner_polygon); 


    unsigned int number_points = 0; 
    for (unsigned int i = 0; i < outer_polygon.size(); ++i)
    {
        number_points += outer_polygon[i].size(); 
    }
    for (unsigned int i = 0; i < inner_polygon.size(); ++i)
    {
        number_points += inner_polygon[i].size(); 
    }

    std::vector<std::vector<double>> coords(number_points, std::vector<double>(2,0)); 
    
    unsigned int index = 0; 
    for (unsigned int i = 0; i < outer_polygon.size(); ++i)
    {
        for (unsigned int j = 0; j < outer_polygon[i].size(); ++j)
        {
            coords[index  ][0] = outer_polygon[i][j][0]; // x-coordinate
            coords[index++][1] = outer_polygon[i][j][1]; // y-coordinate
        }
    }
    for (unsigned int i = 0; i < inner_polygon.size(); ++i)
    {
        for (unsigned int j = 0; j < inner_polygon[i].size(); ++j)
        {
            coords[index  ][0] = inner_polygon[i][j][0]; // x-coordinate
            coords[index++][1] = inner_polygon[i][j][1]; // y-coordinate
        }
    }

    return coords;     
}

std::vector<std::vector<double>> EmbeddedIgaModeler::PrintParametricTriangulation()
{
    const auto face = m_brep_model_vector[0].GetFaceVector()[0];
    std::vector<std::vector<array_1d<double, 2>>> outer_polygon;
    std::vector<std::vector<array_1d<double, 2>>> inner_polygon;
    
    EmbeddedIgaTessellation::CreateTessellation(
        face, outer_polygon, inner_polygon); 

    std::vector<Matrix> triangulation;

    EmbeddedIgaTriangulation embedded_triangulation; 
    embedded_triangulation.CreateTriangulation(
        outer_polygon, inner_polygon, triangulation);

    std::vector<std::vector<double>> coords(triangulation.size() * 3, std::vector<double>(2,0)); 
        
    unsigned int point_index = 0; 
    for (unsigned int tri_i = 0; tri_i < triangulation.size(); ++tri_i)
    {
        for (unsigned int j = 0; j < 3; ++j)    
        {
            for (unsigned int k = 0; k < 2; ++k)    coords[point_index][k] = triangulation[tri_i](j,k); 
            point_index += 1; 
        }
    }

    return coords; 
}

std::vector<std::vector<double>> EmbeddedIgaModeler::PrintMappedPoints()
{
    const auto face = m_brep_model_vector[0].GetFaceVector()[0];
    std::vector<std::vector<array_1d<double, 2>>> outer_polygon_uv;
    std::vector<std::vector<array_1d<double, 2>>> inner_polygon_uv;
    
    EmbeddedIgaTessellation::CreateTessellation(
        face, outer_polygon_uv, inner_polygon_uv); 

    std::vector<Matrix> triangulation_uv;

    EmbeddedIgaTriangulation embedded_triangulation; 
    embedded_triangulation.CreateTriangulation(
        outer_polygon_uv, inner_polygon_uv, triangulation_uv);


    std::vector<Matrix> triangulation_xyz;
    MapCartesianSpace(
        face, triangulation_uv, triangulation_xyz);

    std::vector<std::vector<double>> coords(triangulation_xyz.size() * 3, std::vector<double>(3,0)); 
        
    unsigned int point_index = 0; 
    for (unsigned int tri_i = 0; tri_i < triangulation_xyz.size(); ++tri_i)
    {
        for (unsigned int j = 0; j < 3; ++j)    
        {
            for (unsigned int k = 0; k < 3; ++k)    coords[point_index][k] = triangulation_xyz[tri_i](j,k); 
            point_index += 1; 
        }
    }
    return coords; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
    m_model_part(rModelPart)
{
}
} // namespace Kratos.