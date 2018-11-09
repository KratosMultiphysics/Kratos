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
            if (m_brep_edges[i].Id() == brep_id)
            {
                return success;
                m_brep_edges[i].GetGeometryIntegration(rModelPart, rType, rName, rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);
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
}  // namespace Kratos.

