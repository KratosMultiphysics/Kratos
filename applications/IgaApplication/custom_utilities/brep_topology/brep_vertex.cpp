//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// Project includes
#include "brep_vertex.h"

namespace Kratos
{
    const BrepVertex::VertexTopology BrepVertex::GetVertexTopology(
        const int& rTopologyIndex) const
    {
        KRATOS_ERROR_IF(rTopologyIndex >= m_BrepVertexTopologyVector.size())
            << "Number of topology references smaller than selected index!" << std::endl;

        return m_BrepVertexTopologyVector[rTopologyIndex];
    }
}  // namespace Kratos.

