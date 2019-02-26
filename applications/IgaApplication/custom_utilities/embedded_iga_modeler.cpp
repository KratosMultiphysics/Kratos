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
void EmbeddedIgaModeler::CreateElements3D(ModelPart& rSkinModelPart)
{
    for(unsigned int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
    {
        for (unsigned int face_i = 0; face_i < m_brep_model_vector[brep_i].GetFaceVector().size(); ++face_i)
        {
            const auto face = m_brep_model_vector[brep_i].GetFaceVector()[face_i];

            std::vector<std::vector<array_1d<double, 2>>> outer_polygon_uv;
            std::vector<std::vector<array_1d<double, 2>>> inner_polygon_uv;
            EmbeddedIgaTessellation::CreateTessellation(
                face, outer_polygon_uv, inner_polygon_uv); 

            EmbeddedIgaTriangulation embedded_triangulation; 

            std::vector<Matrix> triangulation_xyz;
            embedded_triangulation.CreateTriangulation(
                face, outer_polygon_uv, inner_polygon_uv, triangulation_xyz);


            unsigned int node_id = 0;
            for (unsigned int tri_i = 0; tri_i < triangulation_xyz.size(); ++tri_i)
            {
                for (unsigned int point_i = 0; point_i < 3; ++point_i)
                {
                    rSkinModelPart.CreateNewNode(node_id, 
                        triangulation_xyz[tri_i](point_i,0), 
                        triangulation_xyz[tri_i](point_i,1), 
                        triangulation_xyz[tri_i](point_i,2));
                    node_id += 1; 
                }
            }

            Properties::Pointer prop = rSkinModelPart.pGetProperties(0);

            node_id = 0;
            // create elements in skin_model_part
            for (unsigned int element_i = 0; element_i < triangulation_xyz.size(); ++element_i)
            {
                rSkinModelPart.CreateNewElement("Element3D3N", element_i, {{node_id, node_id + 1, node_id + 2}}, prop);
                node_id += 3;
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
        face, outer_polygon, inner_polygon, triangulation);

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

    EmbeddedIgaTriangulation embedded_triangulation; 
    std::vector<Matrix> triangulation_uv;
    embedded_triangulation.CreateTriangulation(
        face, outer_polygon_uv, inner_polygon_uv, triangulation_uv);
        
    std::vector<Matrix> triangulation_xyz;
    EmbeddedIgaMapper::MapCartesianSpace(
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