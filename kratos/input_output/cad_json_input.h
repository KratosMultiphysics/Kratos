//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_CAD_JSON_INPUT_INCLUDED )
#define  KRATOS_CAD_JSON_INPUT_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Geometries
#include "geometries/coupling_geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "geometries/brep_surface.h"
#include "geometries/brep_curve_on_surface.h"

namespace Kratos
{

///@name Kratos Classes
///@{
/// Input for CAD-files.
/** Gives IO capabilities for Nurbs based Brep models in the JSON format defined in
https://amses-journal.springeropen.com/articles/10.1186/s40323-018-0109-4. */
template<class TNodeType = Node<3>, class TEmbeddedNodeType = Point>
class CadJsonInput : public IO
{
    public:

        ///@}
        ///@name Type Definitions
        ///@{

        /// Pointer definition of CadJsonInput
        KRATOS_CLASS_POINTER_DEFINITION(CadJsonInput);

        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Geometry<TNodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;

        typedef PointerVector<TNodeType> ContainerNodeType;
        typedef PointerVector<TEmbeddedNodeType> ContainerEmbeddedNodeType;

        typedef CouplingGeometry<TNodeType> CouplingGeometryType;

        typedef NurbsSurfaceGeometry<3, ContainerNodeType> NurbsSurfaceType;
        typedef NurbsCurveGeometry<2, ContainerEmbeddedNodeType> NurbsTrimmingCurveType;

        typedef typename NurbsSurfaceType::Pointer NurbsSurfacePointerType;
        typedef typename NurbsTrimmingCurveType::Pointer NurbsTrimmingCurvePointerType;

        typedef BrepSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepSurfaceType;
        typedef BrepCurveOnSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;

        typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceLoopType;
        typedef DenseVector<DenseVector<typename BrepCurveOnSurfaceType::Pointer>> BrepCurveOnSurfaceLoopArrayType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor with path to input file.
        CadJsonInput(
            const std::string& rDataFileName,
            SizeType EchoLevel = 0)
            : mEchoLevel(EchoLevel)
        {
            mCadJsonParameters = ReadParamatersFile(rDataFileName, EchoLevel);
        }

        /// Constructor with KratosParameters.
        CadJsonInput(
            Parameters CadJsonParameters,
            SizeType EchoLevel = 0)
            : mCadJsonParameters(CadJsonParameters)
            , mEchoLevel(EchoLevel)
        {
        }

        /// Destructor.
        ~CadJsonInput() = default;

        ///@}
        ///@name Python exposed Functions
        ///@{

        /// Adds all CAD geometries to the herin provided model_part.
        void ReadModelPart(ModelPart& rModelPart) override;

        ///@}

    private:
        ///@name Static Functions
        ///@{

        /// Allows static access without own memory.
        static void ReadGeometryModelPart(
            const Parameters& rCadJsonParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Read in Brep
        ///@{

        static void ReadBreps(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        static void ReadBrepFaces(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        static void ReadBrepEdges(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Read in Brep Geometries
        ///@{

        static void ReadBrepSurfaces(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        static void ReadBrepSurface(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Read in Surface Trimming
        ///@{

        static BrepCurveOnSurfaceLoopType
            ReadTrimmingCurveVector(
                const Parameters& rParameters,
                typename NurbsSurfaceType::Pointer pNurbsSurface,
                ModelPart& rModelPart,
                SizeType EchoLevel = 0);

        static typename BrepCurveOnSurfaceType::Pointer
            ReadTrimmingCurve(
                const Parameters& rParameters,
                typename NurbsSurfaceType::Pointer pNurbsSurface,
                ModelPart& rModelPart,
                SizeType EchoLevel = 0);

        static std::tuple<BrepCurveOnSurfaceLoopArrayType, BrepCurveOnSurfaceLoopArrayType>
            ReadBoundaryLoops(
                const Parameters& rParameters,
                typename NurbsSurfaceType::Pointer pNurbsSurface,
                ModelPart& rModelPart,
                SizeType EchoLevel = 0);

        ///@}
        ///@name Read in Nurbs Geometries
        ///@{

        static void ReadBrepCurveOnSurfaces(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        static void ReadBrepEdge(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        static void ReadBrepEdgeBrepCurveOnSurface(
            const Parameters & rParameters,
            ModelPart & rModelPart,
            SizeType EchoLevel = 0);

        static void ReadCouplingGeometry(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Read in Nurbs Geometries
        ///@{

        /* @brief read NurbsCurves from the given parameter input.
        *
        * The input needs to be provided in the following shape:
        * {
        *     "is_rational": bool,
        *     "degree": p:int,
        *     "knot_vector": [ double, ... ],
        *     "active_range": [ double, double ],
        *     "control_points": [
        *         [ id: int, [ x, y, z, weight ] ],
        *         ...
        *     ]
        * }
        */
        template<int TWorkingSpaceDimension, class TThisNodeType>
        static typename NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer
            ReadNurbsCurve(
                const Parameters& rParameters,
                ModelPart& rModelPart,
                SizeType EchoLevel = 0)
        {
            bool is_rational = true;
            if (rParameters.Has("is_rational")) {
                is_rational = rParameters["is_rational"].GetBool();
            }
            else {
                KRATOS_INFO_IF("ReadNurbsCurve", (EchoLevel > 4))
                    << "\"is_rational\" is not provided within \"surface\". Thus, it is considered as rational. "
                    << "If this curve is non-rational the computation of the shape functions is less optimized."
                    << std::endl;
            }

            KRATOS_ERROR_IF_NOT(rParameters.Has("knot_vector"))
                << "Missing 'knot_vector' in nurbs curve" << std::endl;
            Vector knot_vector = rParameters["knot_vector"].GetVector();

            KRATOS_ERROR_IF_NOT(rParameters.Has("degree"))
                << "Missing 'degree' in nurbs curve" << std::endl;
            int polynomial_degree = rParameters["degree"].GetInt();

            PointerVector<TThisNodeType> control_points;

            ReadControlPointVector(control_points,
                rParameters["control_points"], rModelPart, EchoLevel);

            if (is_rational)
            {
                Vector control_point_weights = ReadControlPointWeightVector(
                    rParameters["control_points"]);

                return Kratos::make_shared<NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>>(
                    NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>(
                        control_points,
                        polynomial_degree,
                        knot_vector));
            }
            return Kratos::make_shared<NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>>(
                NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>(
                    control_points,
                    polynomial_degree,
                    knot_vector));
        }

        /* @brief read NurbsSurfaces from the given parameter input.
        *
        * The input needs to be provided in the following shape:
        * {
        *     "is_trimmed": bool,
        *     "is_rational" : bool,
        *     "degrees" : [int, int] ,
        *     "knot_vectors" : [
        *         [ double, ... ],
        *             [double, ...]
        *     ],
        *     "control_points" : [
        *         [ id:p, [0, 10, 0, 1] ],
        *         ...
        *     ]
        * }
        */
        template<int TWorkingSpaceDimension, class TThisNodeType>
        static typename NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer
            ReadNurbsSurface(
                const Parameters& rParameters,
                ModelPart& rModelPart,
                SizeType EchoLevel = 0)
        {
            bool is_rational = true;
            if (rParameters.Has("is_rational")) {
                is_rational = rParameters["is_rational"].GetBool();
            }
            else {
                KRATOS_INFO_IF("ReadNurbsSurface", (EchoLevel > 4))
                    << "\"is_rational\" is not provided within \"surface\". Thus, it is considered as rational. "
                    << "If this surface is non-rational the computation of the shape functions is less optimized."
                    << std::endl;
            }

            KRATOS_ERROR_IF_NOT(rParameters.Has("knot_vectors"))
                << "Missing 'knot_vector' in nurbs surface" << std::endl;
            KRATOS_ERROR_IF(rParameters["knot_vectors"].size() != 2)
                << "'knot_vectors' need to be of size two, knot_vector_u and knot_vector_v" << std::endl;
            Vector knot_vector_u = rParameters["knot_vectors"][0].GetVector();
            Vector knot_vector_v = rParameters["knot_vectors"][1].GetVector();

            KRATOS_ERROR_IF_NOT(rParameters.Has("degrees"))
                << "Missing 'degrees' in nurbs surface" << std::endl;
            KRATOS_ERROR_IF(rParameters["degrees"].size() != 2)
                << "'degrees' need to be of size two, p and q" << std::endl;
            int p = rParameters["degrees"][0].GetInt();
            int q = rParameters["degrees"][1].GetInt();

            PointerVector<TThisNodeType> control_points;

            ReadControlPointVector(control_points,
                rParameters["control_points"], rModelPart, EchoLevel);

            if (is_rational)
            {
                Vector control_point_weights = ReadControlPointWeightVector(
                    rParameters["control_points"]);

                return Kratos::make_shared<NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>>(
                    control_points,
                    p,
                    q,
                    knot_vector_u,
                    knot_vector_v,
                    control_point_weights);
            }
            return Kratos::make_shared<NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>>(
                NurbsSurfaceGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>(
                    control_points,
                    p,
                    q,
                    knot_vector_u,
                    knot_vector_v));
        }

        ///@}
        ///@name Read in Control Points
        ///@{

        /// Reads the weights of all control points and provides them in a Vector.
        static Vector ReadControlPointWeightVector(
            const Parameters& rParameters,
            SizeType EchoLevel = 0);

        /// Reads a Node<3>::Pointer-vector of control points.
        static void ReadControlPointVector(
            PointerVector<Node<3>>& rControlPoints,
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        /// Reads a Point::Pointer-vector of control points.
        static void ReadControlPointVector(
            PointerVector<Point>& rControlPoints,
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Read in Nodes/ Points
        ///@{

        /* Reads, and returns a Pointer to Node<3>.
        * Input needs to be a Parameter object:
        * [id, [x, y, z, weight]]
        */
        static Node<3>::Pointer ReadNode(
            const Parameters& rParameters,
            ModelPart& rModelPart,
            SizeType EchoLevel = 0);

        /* Reads, and returns a Pointer to Point.
        * Input needs to be a Parameter object:
        * [[x, y, z, weight]] or [id, [x, y, z, weight]]
        */
        static Point::Pointer ReadPoint(
            const Parameters& rParameters,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Utility functions
        ///@{

        /// Sets the geometry Id with either the 'brep_id' or the 'brep_name'.
        template<class TGeometry>
        static void SetIdOrName(
            const Parameters& rParameters,
            typename TGeometry::Pointer pGeometry)
        {
            if (rParameters.Has("brep_id")) {
                pGeometry->SetId(rParameters["brep_id"].GetInt());
            }
            else if (rParameters.Has("brep_name")) {
                pGeometry->SetId(rParameters["brep_name"].GetString());
            }
        }

        /// Returns the string of either the 'brep_id' or the 'brep_name'. Used for output massages.
        static std::string GetIdOrName(
            const Parameters& rParameters);

        /// Checks if one of the 'brep_id' or the 'brep_name' is provided.
        static bool HasIdOrName(
            const Parameters& rParameters);

        /// Returns the geometry with either the 'brep_id' or the 'brep_name'.
        static typename GeometryType::Pointer GetGeometry(
            const Parameters& rParameters,
            ModelPart& rModelPart);

        /// Reads in a json formatted file and returns its KratosParameters instance.
        static Parameters ReadParamatersFile(
            const std::string& rDataFileName,
            SizeType EchoLevel = 0);

        ///@}
        ///@name Members
        ///@{

        Parameters mCadJsonParameters;
        int mEchoLevel;

        ///@}
    }; // Class CadJsonInput
}  // namespace Kratos.

#endif // KRATOS_CAD_JSON_INPUT_INCLUDED  defined