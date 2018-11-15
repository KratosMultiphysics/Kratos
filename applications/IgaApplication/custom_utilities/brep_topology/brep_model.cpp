//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:    BSD License
//              Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//                   Thomas Oberbichler
//

// Project includes
#include "brep_model.h"

#include "iga_application.h"
#include "iga_application_variables.h"


namespace Kratos
{
    bool BrepModel::GetIntegrationDomain(ModelPart& rModelPart, int& brep_id,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        bool success = false;
        
        for (int i = 0; i < m_brep_edges.size(); ++i)
        {
            std::cout << "ID: " << m_brep_edges[i].Id() << ", brep_id: " << brep_id << std::endl;
            if (m_brep_edges[i].Id() == brep_id)
            {
                m_brep_edges[i].GetGeometryIntegration(rModelPart, rType, rName, rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);
                return true;
            }
        }

        for (int i = 0; i < m_brep_faces.size(); ++i)
        {
            std::cout << "ID: " << m_brep_faces[i].Id() << ", brep_id: " << brep_id << std::endl;
            if (m_brep_faces[i].Id() == brep_id)
            {
                m_brep_faces[i].GetGeometryIntegration(rModelPart, rType, rName, rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);
                return true;
            }
        }

        return success;
    }

    // --------------------------------------------------------------------------
    std::vector<BrepFace>& BrepModel::GetFaceVector()
    {
        return m_brep_faces;
    }
    std::vector<BrepEdge>& BrepModel::GetEdgeVector()
    {
        return m_brep_edges;
    }
    std::vector<BrepVertex>& BrepModel::GetVertexVector()
    {
        return m_brep_vertices;
    }

    BrepModel::BrepModel(
        int& brep_id,
        double& model_tolerance,
        std::vector<BrepFace>& faces,
        std::vector<BrepEdge>& edges,
        std::vector<BrepVertex>& vertices)
        : m_model_tolerance(model_tolerance),
        m_brep_faces(faces),
        m_brep_edges(edges),
        m_brep_vertices(vertices),
        IndexedObject(brep_id),
        Flags()
    {
        std::cout << "m_brep_edges size: " << m_brep_edges.size() << std::endl;
        std::cout << "m_brep_faces size: " << m_brep_faces.size() << std::endl;
    }
}  // namespace Kratos.

