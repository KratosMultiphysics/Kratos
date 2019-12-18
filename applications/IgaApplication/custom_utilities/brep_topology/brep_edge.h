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
    /// Edge.
    /** Edge can deal as global independent Curve, as Brep and as Coupling Brep.
    *   The topology of each edge is given by the TopologyEdge
    */
    class BrepEdge : public IndexedObject, public Flags
    {
    public:
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

        /* Used to embedd additional parts to the structure.*/
        struct EmbeddedPoint {
            int trim_index;
            double local_parameter;

            EmbeddedPoint(const int& rTrimIndex, const double& rLocalParameter)
            {
                trim_index = rTrimIndex;
                local_parameter = rLocalParameter;
            }
        };

        /* Used to trim the curve.*/
        struct TrimmingRange
        {
            int trim_index;
            ANurbs::Interval<double> range;

            TrimmingRange(const int& rTrimIndex, const Vector& rRange)
            {
                trim_index = rTrimIndex;
                range = ANurbs::Interval<double>(rRange[0], rRange[1]);
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

        /* If edge has more then one brep topology object, then the edge can deal as
        *  a coupling edge.
        */
        bool IsCouplingEdge();

        void GetGeometryNodes(
            ModelPart& rModelPart,
            const int& rT) const;

        void GetGeometryVariationNodes(
            ModelPart& rModelPart,
            const int& rT) const;

        void GetIntegrationGeometry(
            ModelPart& rModelPart,
            const std::string& rType,
            const std::string& rName,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables) const;

        void GetIntegrationBrep(
            ModelPart& rModelPart,
            const int& trim_index,
            const std::string& rType,
            const std::string& rName,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables) const;

        
        const EdgeTopology GetEdgeTopology(
            const int rTopologyIndex) const;

        const Kratos::shared_ptr<NodeCurveGeometry3D> GetCurve3d() const;

        const std::vector<EdgeTopology>& GetBrepEdgeTopologyVector() const; 

        const int GetNumberOfEdgeTopologies() const;

        ///Constructor
        BrepEdge(
            const int rBrepId,
            std::vector<EdgeTopology>& rBrepEdgeTopologyVector,
            std::vector<TrimmingRange>& rTrimmingRangeVector,
            std::vector<EmbeddedPoint>& rEmbeddedPoints,
            const int rDegree,
            Vector& rKnotVector,
            Vector& rActiveRange,
            std::vector<int>& rControlPointIds,
            ModelPart& rModelPart);

        /// Destructor.
        virtual ~BrepEdge() {};

    private:
        // topology parameter
        std::vector<EdgeTopology>              mBrepEdgeTopologyVector;
        std::vector<TrimmingRange>             mTrimmingRangeVector;
        std::vector<EmbeddedPoint>             mEmbeddedPoints;

        //anurbs variables
        std::shared_ptr<NodeCurveGeometry3D>   mNodeCurveGeometry3D;

        //3d curve parameters
        Vector                        mActiveRange;
        ModelPart&                    mModelPart;

    }; // Class BrepEdge
}  // namespace Kratos.

#endif // KRATOS_BREP_EDGE_H_INCLUDED defined