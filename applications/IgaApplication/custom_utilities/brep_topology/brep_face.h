#if !defined(KRATOS_BREP_FACE_H_INCLUDED )
#define  KRATOS_BREP_FACE_H_INCLUDED

#pragma once
// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "iga_application_variables.h"
#include "custom_utilities/node_surface_geometry_3d.h"
#include "brep_trimming_curve.h"
#include "brep_boundary_loop.h"

// Kratos includes
#include "includes/model_part.h"

namespace Kratos
{
    ///@name Type Definitions
    ///@{
    typedef std::vector<int> IntVector;
    ///@}
    ///@name Kratos Classes
    ///@{
    /// Short class definition.
    /** Provides topology functions for Face objects.
    */
    class BrepFace : public IndexedObject, public Flags
    {
    public:
        ///@name Type Definitions
        ///@{
        /* Geometry Refinement Parameters are used to pass the full refinement information
        *  needed for face-refinement.
        */
        struct GeometryRefinementParametersSurface
        {
            /// knot variables
            Vector knot_insertions_u;
            Vector knot_insertions_v;
            int multiply_knots_u;
            int multiply_knots_v;
            double max_element_size_u;
            double max_element_size_v;
            /// degree variables
            int order_elevation_p;
            int order_elevation_q;
            int min_order_p;
            int min_order_q;

            /* Constructor */
            GeometryRefinementParametersSurface() {
                knot_insertions_u = ZeroVector(0);
                knot_insertions_v = ZeroVector(0);

                multiply_knots_u = 0;
                multiply_knots_v = 0;

                max_element_size_u = 0.0;
                max_element_size_v = 0.0;

                order_elevation_p = 0;
                order_elevation_q = 0;

                min_order_p = 0;
                min_order_p = 0;
            }
        };

        struct EmbeddedPoint {
            int trim_index;
            Vector local_coordinates;

            EmbeddedPoint(const int& rTrimIndex, const Vector& rLocalCoordinates)
            {
                trim_index = rTrimIndex;
                local_coordinates = rLocalCoordinates;
            }
        };

        /// Pointer definition of KratosNurbsBrepApplication
        KRATOS_CLASS_POINTER_DEFINITION(BrepFace);
        ///@}
        ///@name Life Cycle
        ///@{

        void GetGeometryIntegration(
            ModelPart& rModelPart,
            const std::string& rType,
            const std::string& rName,
            const int& rPropertiesId,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables);

        void GetGeometryIntegrationTrimmed(
            ModelPart& rModelPart,
            const std::string& rType,
            const std::string& rName,
            const int& rPropertiesId,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables);

        void CreateIntegrationElementsConditions(
            std::vector<ANurbs::IntegrationPoint2<double>> rIntegrationPoints,
            ModelPart& rModelPart,
            const std::string& rType,
            const std::string& rName,
            const int& rPropertiesId,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables) const;

        void GetIntegrationBrepEdge(
            ModelPart& rModelPart,
            const int& trim_index,
            const std::string& rType,
            const std::string& rName,
            const int& rPropertiesId,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables);

        void EvaluatePoint(
            const double& rU,
            const double& rV,
            Element::GeometryType::PointsArrayType& rNonZeroControlPoints,
            Vector& rShapeFunction,
            Matrix& rShapeFunctionDerivative,
            Matrix& rShapeFunctionSecondDerivative
            ) const;


        const Kratos::shared_ptr<Curve<2>> GetTrimCurve(
            const int& trim_index) const;

        /*Returns the member NodeSurfaceGeometry3D object of the brep face.
        This object gives a NURBS representation of the surface of the face.*/
        const Kratos::shared_ptr<NodeSurfaceGeometry3D> GetSurface() const;


        /// Constructor.
        BrepFace(
            int rBrepId,
            bool rIsTrimmed,
            bool rIsRational,
            std::vector<BrepBoundaryLoop>& rTrimmingLoops,
            std::vector<BrepBoundaryLoop>& rEmbeddedLoops,
            std::vector<EmbeddedPoint>& rEmbeddedPoints,
            Vector& rKnotVectorU,
            Vector& rKnotVectorV,
            int& rP,
            int& rQ,
            IntVector& rControlPointIds,
            ModelPart& rModelPart);

        /// Destructor.
        virtual ~BrepFace() {};

        ///@}
    protected:

    private:

        ///@name Private methods
        ///@{

        ///@}
        ///@name Member Variables
        ///@{
        bool m_is_trimmed;
        bool m_is_rational;
        std::vector<BrepBoundaryLoop> m_trimming_loops;
        std::vector<BrepBoundaryLoop> m_embedded_loops;
        std::vector<EmbeddedPoint> m_embedded_points;

        //3d surface parameters
        IntVector& m_control_points_ids;
        ModelPart& m_model_part;

        //anurbs variables
        std::shared_ptr<NodeSurfaceGeometry3D> m_node_surface_geometry_3d;
        ///@}

    }; // Class BrepFace

} // namespace Kratos.

#endif // KRATOS_BREP_FACE_H_INCLUDED  defined