//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Dagmawi Bekel
//

// System includes

// External includes

// Project includes
#include "embedded_iga_modeler.h"

namespace Kratos
{
void EmbeddedIgaModeler::CreateTessellation(std::shared_ptr<NodeCurveGeometry3D> nodes_vector)
{
    for (int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
    {
        for (int edge_i = 0; edge_i < m_brep_model_vector[brep_i].GetEdgeVector().size(); ++edge_i)
        {
            const auto node = m_brep_model_vector[0].GetEdgeVector()[0].GetCurve3d()->GetNode(0);
            
            nodes_vector->SetNode(edge_i,node); 
        }
    }
}

void EmbeddedIgaModeler::CreateElements()
{
    auto nodes_vector = CreateTessellation();


}

EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart &rModelPart)
    : NurbsBrepModeler::NurbsBrepModeler(rModelPart),
      m_model_part(rModelPart)
{
}
} // namespace Kratos.
