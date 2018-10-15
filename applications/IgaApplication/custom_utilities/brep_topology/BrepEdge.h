//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//                   Thomas Oberbichler
//


#if !defined(KRATOS_BREP_EDGE_H_INCLUDED )
#define  KRATOS_BREP_EDGE_H_INCLUDED


// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

#include "../node_curve_geometry_3d.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /// Edge in global space.
    /** Detail class definition.
    */
    class BrepEdge : public IndexedObject, public Flags
    {
    public:
        /* Used to separate the curve into trimmed ranges. */
        struct TrimmingRange
        {
            int trim_index;
            Vector range;

            TrimmingRange(const int& rTrimIndex, const Vector& rRange)
            {
                trim_index = rTrimIndex;
                range = rRange;
            }
        };

        /* Used to descibe the topology of edges. */
        struct EdgeTopology
        {
            int face_id;
            int trim_index;
            bool relative_direction;

            Topology(const int& rFaceId, const int& rTrimIndex, const bool& rRelativeDirection)
            {
                face_id = rFaceId;
                trim_index = rTrimIndex;
                relative_direction = rRelativeDirection;
            }
        };

        ///@name Life Cycle 
        ///@{ 

        bool IsCouplingEdge();

        ///Constructor
        BrepEdge::BrepEdge(unsigned int edge_id,
            std::vector<Topology>& brep_edge_topology_vector,
            std::vector<TrimmingRange>& trimming_range_vector,
            unsigned int& degree,
            Vector& knot_vector,
            Vector& active_range,
            std::vector<int>& control_point_ids,
            Kratos::shared_ptr<ModelPart> model_part);

        /// Destructor.
        virtual ~BrepEdge() {};

        ///@} 
    protected:

    private:

        ///@name Member Variables
        ///@{ 
        // topology parameter
        std::vector<EdgeTopology>  m_brep_edge_topology_vector;
        std::vector<TrimmingRange> m_trimming_range_vector;

        NodeCurveGeometry3D        m_curve_geometry;

        //3d curve parameter
        unsigned int                  m_degree;
        Vector                        m_knot_vector;
        Vector                        m_active_range;
        std::vector<int>              m_control_point_ids;
        Kratos::shared_ptr<ModelPart> mp_model_part;
        ///@}    
    }; // Class BrepEdge 

}  // namespace Kratos.

#endif // KRATOS_BREP_EDGE_H_INCLUDED  defined