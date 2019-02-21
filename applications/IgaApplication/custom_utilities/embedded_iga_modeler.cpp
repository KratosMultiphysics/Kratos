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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
    m_model_part(rModelPart)
{
}
} // namespace Kratos.