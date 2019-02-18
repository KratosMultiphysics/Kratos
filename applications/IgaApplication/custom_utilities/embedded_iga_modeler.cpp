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
std::vector<std::vector<double>> EmbeddedIgaModeler::Test()
{
    const auto face = m_brep_model_vector[0].GetFaceVector()[0];

    EmbeddedIgaTessellation tessellation;
    auto tes = tessellation.CreateTessellation(face);




    // std::vector<std::vector<array_1d<double,2>>> poly; 
    // std::vector<array_1d<double,2>> loop(3); 
    
    // loop[0][0] = 1; 
    // loop[0][1] = 2; 
    // loop[1][0] = 3; 
    // loop[1][1] = 4; 
    // loop[2][0] = 5; 
    // loop[2][1] = 6; 

    // poly.push_back(loop);

    // KRATOS_WATCH(loop.size())
    // KRATOS_WATCH(poly.size())
    // KRATOS_WATCH(poly[0].size())
    // KRATOS_WATCH(poly)

    
    return tes;
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