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

// System includes

// Project includes
#include "includes/model_part.h"
#include "iga_application_variables.h"
#include "custom_utilities/node_curve_geometry_3d.h"


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

        //using Kratos::NurbsBrepModeler;

        /// Pointer definition of KratosNurbsBrepApplication
        KRATOS_CLASS_POINTER_DEFINITION(BrepEdge);

        /* Geometry Refinement Parameters are used to pass the full refinement information
        *  needed for face-refinement.
        */
        struct GeometryRefinementParametersCurve
        {
            /// knot variables
            Vector knot_insertions_u;
            int multiply_knots_u;
            double max_element_size_u;
            /// degree variables
            int order_elevation_p;
            int min_order_p;

            /* Constructor */
            GeometryRefinementParametersCurve() {
                knot_insertions_u = ZeroVector(0);
                multiply_knots_u = 0;
                max_element_size_u = 0.0;
                order_elevation_p = 0;
                min_order_p = 0;
            }
        };

        struct EmbeddedPoint {
            int trim_index;
            double local_parameter;

            EmbeddedPoint(const int& rTrimIndex, const double& rLocalParameter)
            {
                trim_index = rTrimIndex;
                local_parameter = rLocalParameter;
            }
        };

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
            int brep_id;
            int trim_index;
            bool relative_direction;

            EdgeTopology(const int& rBrepId, const int& rTrimIndex, const bool& rRelativeDirection)
            {
                brep_id = rBrepId;
                trim_index = rTrimIndex;
                relative_direction = rRelativeDirection;
            }
        };

        ///@name Life Cycle
        ///@{

        bool IsCouplingEdge();

       void GetGeometryIntegration(ModelPart& rModelPart,
           const std::string& rType,
           const std::string& rName,
           const int& rPropertiesId,
           const int& rShapeFunctionDerivativesOrder,
           std::vector<std::string> rVariables);

        ///Constructor
        BrepEdge(
            int rBrepId,
            std::vector<EdgeTopology>& rBrepEdgeTopologyVector,
            std::vector<TrimmingRange>& rTrimmingRangeVector,
            std::vector<EmbeddedPoint>& rEmbeddedPoints,
            int& rDegree,
            Vector& rKnotVector,
            Vector& rActiveRange,
            std::vector<int>& rControlPointIds,
            ModelPart& rModelPart);

        /// Destructor.
        virtual ~BrepEdge() {};

        ///@}
    protected:

    private:

        ///@name Member Variables
        ///@{
        // topology parameter
        std::vector<EdgeTopology>              m_brep_edge_topology_vector;
        std::vector<TrimmingRange>             m_trimming_range_vector;
        std::vector<EmbeddedPoint>             m_embedded_points;

        //anurbs variables
        std::shared_ptr<NodeCurveGeometry3D>   m_node_curve_geometry_3d;

        //3d curve parameters
        Vector                        m_active_range;
        std::vector<int>              m_control_point_ids;
        ModelPart&                    m_model_part;
        ///@}
    }; // Class BrepEdge

}  // namespace Kratos.

#endif // KRATOS_BREP_EDGE_H_INCLUDED defined