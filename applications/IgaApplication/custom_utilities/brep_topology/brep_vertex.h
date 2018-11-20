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
//


#if !defined(KRATOS_BREP_VERTEX_H_INCLUDED )
#define  KRATOS_BREP_VERTEX_H_INCLUDED


// Project includes
#include "iga_application_variables.h"


namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /// Edge in global space.
    /** Vertex in spatial coordinates.
    */
    class BrepVertex : public IndexedObject, public Flags
    {
    public:
        /* Descibes the topology properties of vertices. */
        struct VertexTopology
        {
            int brep_id;
            int trim_index;

            VertexTopology(const int& rBrepId, const int& rTrimIndex)
            {
                brep_id    = rBrepId;
                trim_index = rTrimIndex;
            }
        };

        ///@name Life Cycle
        ///@{
        VertexTopology GetVertexInformation(const int& rTopologyIndex);

        const BrepVertex::VertexTopology BrepVertex::GetVertexTopology(
            const int& rTopologyIndex) const;

        /// Constructor.
        BrepVertex::BrepVertex(unsigned int& rVertexId,
            std::vector<VertexTopology>& rBrepVertexTopologyVector,
            int rControlPointId,
            Vector& rPoint)
            : m_brep_vertex_topology_vector(rBrepVertexTopologyVector),
              m_control_point_id(rControlPointId),
              m_point(rPoint),
              IndexedObject(rVertexId),
              Flags()
        {};

        /// Destructor.
        virtual ~BrepVertex() {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        // topology
        std::vector<VertexTopology> m_brep_vertex_topology_vector;

        //3d point
        // TODO include KRATOS node as spatial coordinates, this would allow to neglect
        // following attributes.
        int    m_control_point_id;
        Vector m_point;
        ///@}
    }; // Class BrepVertex

}  // namespace Kratos.

#endif // KRATOS_BREP_VERTEX_H_INCLUDED  defined