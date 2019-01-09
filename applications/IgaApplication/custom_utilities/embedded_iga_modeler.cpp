//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Thomas Oberbichler
//

// System includes

// External includes

// Project includes
#include "embedded_iga_modeler.h"


namespace Kratos
{
    Matrix EmbeddedIgaModeler::CreateTessellation()
    {
        Matrix Nodes = ZeroMatrix(m_brep_model_vector[0].GetEdgeVector().size(), 6);
        for (int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
        {
            
            for (int edge_i = 0; edge_i < m_brep_model_vector[brep_i].GetEdgeVector().size(); ++edge_i)
            {
                const auto start_node = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(0)->Coordinates();
                const auto end_node = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(1)->Coordinates();
                
                
                Nodes(edge_i,0) = start_node[0];
                Nodes(edge_i,1) = start_node[1];
                Nodes(edge_i,2) = start_node[2];
                Nodes(edge_i,3) = end_node[0];
                Nodes(edge_i,4) = end_node[1];
                Nodes(edge_i,5) = end_node[2];
                

            }
        }
        return Nodes; 
    }



    EmbeddedIgaModeler::EmbeddedIgaModeler(ModelPart& rModelPart)
        :  NurbsBrepModeler::NurbsBrepModeler(rModelPart),
         m_model_part(rModelPart)
    {
    }
}  // namespace Kratos.


