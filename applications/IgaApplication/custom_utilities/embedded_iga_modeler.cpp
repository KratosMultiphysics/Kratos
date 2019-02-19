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