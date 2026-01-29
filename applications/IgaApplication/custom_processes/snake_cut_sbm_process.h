// //    |  /           |
// //    ' /   __| _` | __|  _ \   __|
// //    . \  |   (   | |   (   |\__ `
// //   _|\_\_|  \__,_|\__|\___/ ____/
// //                   Multi-Physics
// //
// //  License:         BSD License
// //                   Kratos default license: kratos/license.txt
// //
// //  Main authors:    Nicolo' Antonelli
// //                   Andrea Gorgi
// //

// #pragma once

// // System includes
// #include <utility>
// #include <vector>

// // Project includes
// #include "containers/model.h"
// #include "includes/model_part.h"
// #include "spatial_containers/bins_dynamic.h"
// #include "processes/process.h"
// #include "geometries/nurbs_curve_geometry.h"
// #include "snake_sbm_process.h"
// #include "custom_utilities/create_breps_sbm_utilities.h"
// #include "includes/global_pointer_variables.h"

// namespace Kratos
// {
// ///@name Kratos Classes
// ///@{

// /**
//  * @class SnakeCutSbmProcess
//  * @brief Process class for implementing Snake-based Surrogate Boundary Method (SBM).
//  * This class provides various functions to create and manipulate surrogate boundaries
//  * for Iga models in Kratos.
//  */
// class KRATOS_API(IGA_APPLICATION) SnakeCutSbmProcess
//     : public SnakeSbmProcess
// {

// public:

//     typedef Node NodeType;
//     typedef PointerVector<Node> ContainerNodeType;
//     typedef PointerVector<Point> ContainerEmbeddedNodeType;
//     typedef BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;
//     typedef typename ModelPart::ElementsContainerType ElementsContainerType;

//     typedef Geometry<NodeType> GeometryType;
//     typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
//     typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;
//     typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
//     typedef typename Properties::Pointer PropertiesPointerType;
//     typedef BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveType;
    
//     using NurbsCurveGeometryType = NurbsCurveGeometry<3, PointerVector<Node>>;
//     typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceType;
//     using ProjectionSegment = std::pair<CoordinatesArrayType, CoordinatesArrayType>;
//     using NodePointerVector = GlobalPointersVector<NodeType>;

//     ///@name Type Definitions
//     ///@{

//     /// Pointer definition of SnakeCutSbmProcess
//     KRATOS_CLASS_POINTER_DEFINITION(SnakeCutSbmProcess);

//     ///@name Life Cycle
//     ///@{

//     /// Constructor
//     SnakeCutSbmProcess(
//         Model& rModel,
//         Parameters ThisParameters);

//     /// Destructor.
//     ~SnakeCutSbmProcess() = default;

//     ///@}
//     ///@name Operations
//     ///@{

//     void Execute() override
//     {   
//         CreateSbmExtendedGeometries();
//     };
    
//     void ExecuteInitialize() override
//     {
//         SnakeSbmProcess::CreateTheSnakeCoordinates();
//     };

//     void ExecuteInitializeSolutionStep() override
//     {};

//     void ExecuteBeforeSolutionLoop() override
//     {};

//     const Parameters GetDefaultParameters() const override
//     {
//         const Parameters default_parameters = Parameters(R"(
//         {
//             "model_part_name" : "",
//             "skin_model_part_inner_initial_name" : "SkinModelPartInnerInitial",
//             "skin_model_part_outer_initial_name" : "SkinModelPartOuterInitial",
//             "skin_model_part_name" : "SkinModelPart",
//             "echo_level" : 0,
//             "lambda_inner" : 0.5,
//             "lambda_outer" : 0.5,
//             "number_of_inner_loops": 0,
//             "number_internal_divisions": 1,
//             "cut_element_name": "",
//             "cut_interface_condition_name": ""
//         })" );

//         return default_parameters;
//     }


// struct BinSearchParameters
// {
//     DynamicBins&                                   TestBins;        //  reference!
//     SizeType                                       NumberOfResults;
//     ModelPart::NodesContainerType::ContainerType   Results;
//     std::vector<double>                            ListOfDistances;
//     double                                         SearchRadius;

//     BinSearchParameters(
//         DynamicBins&                               bins,            // take by ref
//         SizeType                                   n_results,
//         ModelPart::NodesContainerType::ContainerType results,
//         std::vector<double>                        distances,
//         double                                     radius)
//         : TestBins(bins)
//         , NumberOfResults(n_results)
//         , Results(results)
//         , ListOfDistances(distances)
//         , SearchRadius(radius)
//     {}

//     void reset()
//     {
//         Results.clear();
//         ListOfDistances.clear();
//     }
// };


// struct IntegrationParameters
// {
//     SizeType  NumberOfShapeFunctionsDerivatives;
//     IntegrationInfo CurveIntegrationInfo;
//     const Vector&           KnotSpanSizes;          // <- now const
//     IntegrationParameters(
//         SizeType                    n_derivatives,
//         const IntegrationInfo&      info,
//         const Vector&               knot_spans)     // <- const here too
//         : NumberOfShapeFunctionsDerivatives(n_derivatives)
//         , CurveIntegrationInfo(info)
//         , KnotSpanSizes(knot_spans)
//     {}
// };
    
// private:

//     ModelPart* mpCutElementsSubModelPart = nullptr; 
//     ModelPart* mpCutInterfaceSubModelPart = nullptr; 
//     std::string mCutElementName;
//     std::string mCutInterfaceConditionName;
//     SizeType mInternalDivision;
//     SizeType mCutApproximationOrder;

//     void CreateSbmExtendedGeometries();

//     /**
//      * @brief 
//      * 
//      * @tparam TIsInnerLoop 
//      * @param rSkinSubModelPart 
//      * @param rSurrogateSubModelPart 
//      */
//     template <bool TIsInnerLoop>
//     void CreateSbmExtendedGeometries(
//         const ModelPart& rSkinSubModelPart,
//         const ModelPart& rSurrogateSubModelPart);

    
//     template <bool TIsInnerLoop>
//     void CreateCutAndSkinQuadraturePoints(
//         IntegrationParameters& rIntegrationParameters,
//         BinSearchParameters& rBinSearchParameters,
//         NurbsSurfaceType::Pointer& pNurbsSurface,
//         const Node::Pointer& pSurrogateNode1, 
//         const Node::Pointer& pSurrogateNode2, 
//         const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
//         ModelPart& rIgaModelPart,
//         const ModelPart& rSkinSubModelPart);

//     template <bool TIsInnerLoop>
//     void SetSurrogateToSkinProjections(
//         const ModelPart& rSurrogateSubModelPart, 
//         const ModelPart& rSkinSubModelPart,
//         BinSearchParameters& rBinSearchParameters, 
//         BinSearchParameters& rBinSearchInterfaceParameters);
    
//     void AssestProjectionsFeasibility(
//         const ModelPart& rSkinSubModelPart,
//         Node::Pointer pSurrogateNode1, 
//         Node::Pointer pSurrogateNode2,
//         BinSearchParameters& rBinSearchInterfaceParameters);

//     void FindClosestTruePointToSurrogateVertexByNurbs();

//     bool ProjectToSkinBoundary(
//         const ModelPart* pSkinModelPart,
//         const CoordinatesArrayType& rPoint,
//         CoordinatesArrayType& rProjectedPoint,
//         CoordinatesArrayType& rProjectedPointLocal,
//         std::vector<array_1d<double, 3>>& rCurveDerivatives,
//         int nInitialGuesses);

    
//     void CreateConditions(
//         typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
//         typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
//         ModelPart& rModelPart,
//         const std::string& rConditionName,
//         SizeType& rIdCounter,
//         PropertiesPointerType pProperties,
//         const Vector KnotSpanSizes,
//         const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const;

//     /// Creates elements from geometries
//     void CreateElements(
//         typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
//         typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
//         ModelPart& rDestinationModelPart,
//         const std::string& rElementName,
//         SizeType& rIdCounter,
//         PropertiesPointerType pProperties,
//         const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const;

//     static typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer CreateBrepCurve(
//         const Node::Pointer pFirstBrepPoint,
//         const Node::Pointer pSecondBrepPoint, 
//         const Vector& rActiveRangeKnotVector)
//         {
//             // Create the data for the trimming curves
//             PointerVector<Node> control_points;
//             control_points.push_back(pFirstBrepPoint);
//             control_points.push_back(pSecondBrepPoint);
//             const int polynomial_degree = 1;
//             Vector knot_vector = ZeroVector(4) ;
//             knot_vector[0] = rActiveRangeKnotVector[0] ;
//             knot_vector[1] = rActiveRangeKnotVector[0] ;
//             knot_vector[2] = rActiveRangeKnotVector[1] ;
//             knot_vector[3] = rActiveRangeKnotVector[1] ;
//             // Create the trimming curves
//             typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_trimming_curve(
//                 new NurbsCurveGeometry<3, PointerVector<Node>>(
//                     control_points,
//                     polynomial_degree,
//                     knot_vector));   
//             return p_trimming_curve;
//         }
    
//     /// Return n×n Gauss points mapped on the curved quadrilateral
//     IntegrationPointsArrayType CreateCoonsPatchGaussPoints(
//         std::size_t                  Order,
//         const GeometryType&         rB0,
//         const GeometryType&         rL0,
//         const GeometryType&         rL1,
//         const GeometryType&         rB1,
//         const array_1d<double,3>&    rP00,
//         const array_1d<double,3>&    rP01,
//         const array_1d<double,3>&    rP10,
//         const array_1d<double,3>&    rP11) const;

//     // --- static utility helpers ----------------------------------

//     /// 1-D Gauss-Legendre nodes & weights mapped to [0,1]
//     static void GaussLegendreOnUnitInterval(
//         std::size_t              Order,
//         std::vector<double>&     rXi,
//         std::vector<double>&     rWeight);

//     /// Global coordinates of a point on a Brep-curve at parameter t
//     static array_1d<double,3> GlobalPoint(
//         const GeometryType& rCurve,
//         double               T);

//     /// Coons patch mapping X(ξ,η)
//     static array_1d<double,3> CoonsPoint(
//         double                     Xi,
//         double                     Eta,
//         const GeometryType&       rB0,
//         const GeometryType&       rL0,
//         const GeometryType&       rL1,
//         const GeometryType&       rB1,
//         const array_1d<double,3>&  rP00,
//         const array_1d<double,3>&  rP01,
//         const array_1d<double,3>&  rP10,
//         const array_1d<double,3>&  rP11);

//     /// Finite-difference derivative of CoonsPoint
//     static array_1d<double,3> CoonsDerivativeFD(
//         double                     Xi,
//         double                     Eta,
//         bool                       WithRespectToXi,
//         const GeometryType&       rB0,
//         const GeometryType&       rL0,
//         const GeometryType&       rL1,
//         const GeometryType&       rB1,
//         const array_1d<double,3>&  rP00,
//         const array_1d<double,3>&  rP01,
//         const array_1d<double,3>&  rP10,
//         const array_1d<double,3>&  rP11,
//         double                     Step = 1.0e-8);


//     /** @brief Build control-points, knot-vector and weights of a quadratic
//      *         NURBS parabola through three *existing* nodes.
//      *
//      *  @param pNode0  start   node  (P0)
//      *  @param pNodeM  middle  node  (Pm) – parametric coord tm (default 0.5)
//      *  @param pNode2  end     node  (P2)
//      *  @param rCtrlPtsPointerVector  OUT pointer-vector with 3 control nodes
//      *                                (new node created for the middle CP)
//      *  @param rKnots               OUT Vector size 6 : [0 0 0 1 1 1]
//      *  @param rWeights             OUT Vector size 3 : [1 1 1]
//      *  @param tm                   parametric coordinate of Pm in (0,1)
//      */
//     void BuildParabolicNurbsData(
//         Node::Pointer           pNode0,
//         Node::Pointer           pNodeM,
//         Node::Pointer           pNode2,
//         PointerVector<Node>&    rCtrlPtsPointerVector,
//         Vector&                 rKnots,
//         Vector&                 rWeights,
//         const double            tm = 0.5);

//     /** Collects skin points between two condition IDs (inclusive), with wrap-around. */
//     template <bool TIsInnerLoop>
//     std::vector<array_1d<double,3>> CollectSkinPointsBetween(
//         const ModelPart& rSkinSubModelPart,
//         IndexType id_cond_1,
//         IndexType id_cond_2) const
//     {
//         std::vector<array_1d<double,3>> pts;
//         pts.reserve(rSkinSubModelPart.NumberOfConditions());

//         const IndexType first_id = rSkinSubModelPart.ConditionsBegin()->Id();
//         const IndexType last_id  = first_id + rSkinSubModelPart.NumberOfConditions() - 1;
//         auto next_id = [&](IndexType id){ return (id < last_id) ? (id + 1) : first_id; };
//         auto previous_id = [&](IndexType id){ return (id > first_id) ? (id - 1) : last_id; };

//         IndexType id = id_cond_1;
        
//         int iter = 0;
//         while (true) {
//             const auto& cond = rSkinSubModelPart.GetCondition(id);
//             const auto& node = cond.GetGeometry().pGetPoint(0);
//             pts.push_back(node->Coordinates());
//             if (id == id_cond_2) break;

//             if constexpr (TIsInnerLoop)  
//             { // inner loop: go against the condition orientation
//                 id = previous_id(id);
//             } else { // outer loop: follow the condition orientation
//                 id = next_id(id);
//             }
//             iter++;
//             KRATOS_ERROR_IF(iter > rSkinSubModelPart.NumberOfConditions())
//                 << "CollectSkinPointsBetween: infinite loop detected between IDs "
//                 << id_cond_1 << " and " << id_cond_2 << ".\n";
//         }

//         return pts;
//     }

//     /** Chord-length parametrization on [0,1]. */
//     std::vector<double> ChordLengthParams01(const std::vector<array_1d<double,3>>& Q) const
//     {
//         const std::size_t n = Q.size();
//         std::vector<double> t(n,0.0);
//         if (n <= 1) return t;

//         double L = 0.0;
//         std::vector<double> d(n,0.0);
//         for (std::size_t i=1; i<n; ++i) {
//             d[i] = norm_2(Q[i]-Q[i-1]);
//             L += d[i];
//         }
//         if (L <= 1e-16) return t;

//         double acc = 0.0;
//         for (std::size_t i=1; i<n; ++i) {
//             acc += d[i];
//             t[i] = acc / L;
//         }
//         t.front() = 0.0;
//         t.back()  = 1.0;
//         return t;
//     }

//     /** Closed-form least squares for the middle control point P1 (P0,P2 fixed). */
//     array_1d<double,3> SolveP1LeastSquares(
//         const array_1d<double,3>& P0,
//         const array_1d<double,3>& P2,
//         const std::vector<array_1d<double,3>>& Q,
//         const std::vector<double>& t) const
//     {
//         array_1d<double,3> num = ZeroVector(3);
//         double denom = 0.0;

//         for (std::size_t i=0; i<Q.size(); ++i) {
//             const double ti = std::min(1.0, std::max(0.0, t[i]));
//             const double a = (1.0 - ti)*(1.0 - ti);
//             const double b = 2.0 * ti * (1.0 - ti);
//             const double c = ti*ti;

//             const array_1d<double,3> di = Q[i] - (a*P0 + c*P2);
//             num   += b * di;
//             denom += b * b;
//         }

//         if (denom <= 1e-20) {
//             // Degenerate: fallback to the midpoint parabola
//             return 0.5*(P0 + P2);
//         }
//         return num / denom;
//     }

//     /** Builds a quadratic (degree-2) polynomial NURBS with control points P0, P1, P2 (weights = 1). */
//     typename NurbsCurveGeometryType::Pointer MakeQuadraticNurbs(
//         const Node::Pointer& pNode0,
//         const array_1d<double,3>& P1,
//         const Node::Pointer& pNode2) const
//     {
//         PointerVector<Node> ctrl_pts;
//         ctrl_pts.reserve(3);

//         // Virtual control node for the middle control point (not stored in a ModelPart)
//         Node::Pointer pNode1 = Node::Pointer(new Node(1, P1));
//         ctrl_pts.push_back(pNode0);
//         ctrl_pts.push_back(pNode1);
//         ctrl_pts.push_back(pNode2);

//         // Quadratic Bézier knot vector
//         Vector knots(6);
//         knots[0]=knots[1]=knots[2]=0.0; knots[3]=knots[4]=knots[5]=1.0;

//         // Polynomial (all weights = 1)
//         Vector weights(3);
//         weights[0]=weights[1]=weights[2]=1.0;

//         const unsigned int degree = 2;
//         return typename NurbsCurveGeometryType::Pointer(
//             new NurbsCurveGeometryType(ctrl_pts, degree, knots, weights));
//     }

//     /** Fits the best quadratic Bézier/NURBS (degree 2) between two skin points using all inner real points. */
//     template <bool TIsInnerLoop>
//     typename NurbsCurveGeometryType::Pointer FitQuadraticBezierBetween(
//         const ModelPart& rSkinSubModelPart,
//         IndexType id_closest_true_condition_1,
//         IndexType id_closest_true_condition_2,
//         const Node::Pointer& p_skin1_brep_point, // P0
//         const Node::Pointer& p_skin2_brep_point  // P2
//     ) const
//     {
//         // 1) Gather real points in-between (inclusive), with wrap-around handling
//         auto Q = CollectSkinPointsBetween<TIsInnerLoop>(
//             rSkinSubModelPart, id_closest_true_condition_1, id_closest_true_condition_2);

//         // Ensure endpoints match the two given skin points
//         if (!Q.empty()) {
//             Q.front() = p_skin1_brep_point->Coordinates();
//             Q.back()  = p_skin2_brep_point->Coordinates();
//         }

//         // 2) Chord-length params on [0,1]
//         auto t = ChordLengthParams01(Q);

//         // 3) Closed-form LS for P1
//         const array_1d<double,3> P0 = p_skin1_brep_point->Coordinates();
//         const array_1d<double,3> P2 = p_skin2_brep_point->Coordinates();
//         const array_1d<double,3> P1 = SolveP1LeastSquares(P0, P2, Q, t);

//         // 4) Build the quadratic NURBS (polynomial)
//         return MakeQuadraticNurbs(p_skin1_brep_point, P1, p_skin2_brep_point);
//     }

//     /** Returns the same quadratic curve but with reversed orientation (t -> 1 - t), without refitting. */
//     typename NurbsCurveGeometryType::Pointer ReverseQuadraticBezier(
//         const typename NurbsCurveGeometryType::Pointer& p_forward) const
//     {
//         KRATOS_ERROR_IF(!p_forward) << "ReverseQuadraticBezier: forward curve is null." << std::endl;
//         const auto polynomial_degree = p_forward->PolynomialDegree(0);
//         KRATOS_ERROR_IF(polynomial_degree != 2)
//             << "ReverseQuadraticBezier: curve degree must be 2." << std::endl;
//         KRATOS_ERROR_IF(p_forward->size() != 3)
//             << "ReverseQuadraticBezier: expected exactly 3 control points." << std::endl;

//         // Access control points via Geometry's Points() container
//         auto& pts = p_forward->Points();
//         Node::Pointer p0 = pts(0);
//         Node::Pointer p1 = pts(1);
//         Node::Pointer p2 = pts(2);

//         // Rebuild with reversed control-point order: [P2, P1, P0]
//         PointerVector<Node> ctrl_rev;
//         ctrl_rev.reserve(3);
//         ctrl_rev.push_back(p2);
//         ctrl_rev.push_back(p1);
//         ctrl_rev.push_back(p0);

//         Vector knots(6);
//         knots[0]=knots[1]=knots[2]=0.0; knots[3]=knots[4]=knots[5]=1.0;

//         Vector weights(3);
//         weights[0]=weights[1]=weights[2]=1.0;

//         const unsigned int degree = 2;
//         return typename NurbsCurveGeometryType::Pointer(
//             new NurbsCurveGeometryType(ctrl_rev, degree, knots, weights));
//     }


//     // --- UV projection ---------------------------------------------------

// /** Project a 3D point to surface local coordinates (u,v). */
// bool ProjectPointToUV(
//     const NurbsSurfaceType& rSurface,
//     const array_1d<double,3>& X,
//     array_1d<double,3>& rUV,   // use [u,v,*], ignore rUV[2]
//     const double tol = 1e-13,
//     const unsigned max_iter = 50) const
// {
//     // rUV.clear(); rUV.resize(3,false); rUV[0]=0.5; rUV[1]=0.5; rUV[2]=0.0; // initial guess
//     // return rSurface.ProjectionPointGlobalToLocalSpace(X, rUV);

//     rUV = X;
//     return true;
// }

// /** Sample a 3D geometry edge at n samples, project to UV, return ordered UV points. */
// std::vector<array_1d<double,3>> SampleEdgeToUV(
//     const GeometryType& rEdge3D,
//     const NurbsSurfaceType& rSurface,
//     const std::size_t n_samples = 9) const
// {
//     std::vector<array_1d<double,3>> uv;
//     uv.reserve(n_samples);

//     const double t0 = 0;//rEdge3D.DomainInterval().GetT0();
//     const double t1 = 1;//rEdge3D.DomainInterval().GetT1();

//     for (std::size_t i=0; i<n_samples; ++i) {
//         const double t = (n_samples==1)? 0.5 : (t0 + (t1-t0)*(double)i/(double)(n_samples-1));
//         CoordinatesArrayType loc(3,0.0); loc[0]=t;
//         array_1d<double,3> X;
//         rEdge3D.GlobalCoordinates(X, loc);

//         array_1d<double,3> uv_i;
//         const bool ok = ProjectPointToUV(rSurface, X, uv_i);
//         KRATOS_ERROR_IF_NOT(ok) << "SampleEdgeToUV: projection to (u,v) failed.\n";
//         uv.push_back(uv_i); // [u,v,*]
//     }
//     return uv;
// }

// /** Collect skin points between two condition IDs, project them to UV (inclusive). */
// template <bool TIsInnerLoop>
// std::vector<array_1d<double,3>> CollectSkinUVBetween(
//     const ModelPart& rSkinSubModelPart,
//     const NurbsSurfaceType& rSurface,
//     IndexType id_node_1,
//     IndexType id_node_2) const
// {
//     std::vector<array_1d<double,3>> uv_pts;
//     uv_pts.reserve(rSkinSubModelPart.NumberOfConditions());

//     const IndexType first_id = rSkinSubModelPart.NodesBegin()->Id();
//     const IndexType last_id  = first_id + rSkinSubModelPart.NumberOfNodes() - 1;
//     auto next_id = [&](IndexType id){ return (id < last_id) ? (id + 1) : first_id; };
//     auto previous_id = [&](IndexType id){ return (id > first_id) ? (id - 1) : last_id; };

//     IndexType id = id_node_1;

//     int iter = 0;
//     while (true) {
//         const auto& p_node = rSkinSubModelPart.pGetNode(id);
//         array_1d<double,3> uv;
//         const bool ok = ProjectPointToUV(rSurface, p_node->Coordinates(), uv);
//         KRATOS_ERROR_IF_NOT(ok) << "CollectSkinUVBetween: projection failed.\n";
//         uv_pts.push_back(uv);
//         if (id == id_node_2) break;

//         if constexpr (TIsInnerLoop)  
//         { // inner loop: go against the condition orientation
//             id = previous_id(id);
//         } else { // outer loop: follow the condition orientation
//             id = next_id(id);
//         }
//         iter++;
//         KRATOS_ERROR_IF(iter > rSkinSubModelPart.NumberOfConditions())
//             << "CollectSkinUVBetween: infinite loop detected between IDs "
//             << id_node_1 << " and " << id_node_2 << ".\n";
//     }
//     return uv_pts;
// }


// // Reuse your alias: NurbsCurveGeometryType = NurbsCurveGeometry<2, PointerVector<Node>>
// // (Working-space dimension = 2 for UV curves)

// /** Chord-length parametrization on [0,1] for UV points. */
// std::vector<double> ChordLengthParams01_UV(const std::vector<array_1d<double,3>>& UV) const
// {
//     const std::size_t n = UV.size();
//     std::vector<double> t(n,0.0);
//     if (n <= 1) return t;

//     double L = 0.0;
//     std::vector<double> d(n,0.0);
//     for (std::size_t i=1; i<n; ++i) {
//         const double du = UV[i][0]-UV[i-1][0];
//         const double dv = UV[i][1]-UV[i-1][1];
//         d[i] = std::sqrt(du*du+dv*dv);
//         L += d[i];
//     }
//     if (L <= 1e-16) return t;

//     double acc = 0.0;
//     for (std::size_t i=1; i<n; ++i) {
//         acc += d[i];
//         t[i] = acc / L;
//     }
//     t.front() = 0.0;
//     t.back()  = 1.0;
//     return t;
// }

// /** Closed-form LS for the middle UV control point (U1,V1), with endpoints fixed. */
// array_1d<double,3> SolveUV_P1_LeastSquares(
//     const array_1d<double,3>& UV0,
//     const array_1d<double,3>& UV2,
//     const std::vector<array_1d<double,3>>& UV,
//     const std::vector<double>& t) const
// {
//     array_1d<double,3> num = ZeroVector(3); // z unused
//     double denom = 0.0;

//     for (std::size_t i=0; i<UV.size(); ++i) {
//         const double ti = std::min(1.0, std::max(0.0, t[i]));
//         const double a = (1.0 - ti)*(1.0 - ti);
//         const double b = 2.0 * ti * (1.0 - ti);
//         const double c = ti*ti;

//         array_1d<double,3> di = UV[i] - (a*UV0 + c*UV2); // z term stays 0
//         num   += b * di;
//         denom += b * b;
//     }
//     if (denom <= 1e-20) { // degenerate
//         array_1d<double,3> mid = 0.5*(UV0 + UV2);
//         mid[2]=0.0; return mid;
//     }
//     array_1d<double,3> P1 = num/denom;
//     P1[2]=0.0; // enforce 2D
//     return P1;
// }

// /** Build a quadratic UV NURBS (degree 2) with three UV control nodes (weights = 1). */
// typename NurbsCurveGeometryType::Pointer MakeQuadraticNurbsUV(
//     const array_1d<double,3>& UV0,
//     const array_1d<double,3>& UV1,
//     const array_1d<double,3>& UV2) const
// {
//     PointerVector<Node> ctrl;
//     ctrl.reserve(3);
//     // Create virtual nodes (Z=0)
//     Node::Pointer n0 = Node::Pointer(new Node(1, array_1d<double,3>({UV0[0],UV0[1],0.0})));
//     Node::Pointer n1 = Node::Pointer(new Node(1, array_1d<double,3>({UV1[0],UV1[1],0.0})));
//     Node::Pointer n2 = Node::Pointer(new Node(1, array_1d<double,3>({UV2[0],UV2[1],0.0})));
//     ctrl.push_back(n0); ctrl.push_back(n1); ctrl.push_back(n2);

//     Vector knots(6);  knots[0]=knots[1]=knots[2]=0.0; knots[3]=knots[4]=knots[5]=1.0;
//     Vector w(3);      w[0]=w[1]=w[2]=1.0;
//     const unsigned degree = 2;
//     return typename NurbsCurveGeometryType::Pointer(
//         new NurbsCurveGeometryType(ctrl, degree, knots, w));
// }

// /** Fit a quadratic UV curve between two skin conditions using all real points in-between. */
// template <bool TIsInnerLoop>
// typename NurbsCurveGeometryType::Pointer FitUV_BetweenSkinConditions(
//     const ModelPart& rSkinSubModelPart,
//     const NurbsSurfaceType& rSurface,
//     IndexType id_cond_1,
//     IndexType id_cond_2) const
// {
//     // UV samples along the true skin between the two conditions
//     auto UV = CollectSkinUVBetween<TIsInnerLoop>(rSkinSubModelPart, rSurface, id_cond_1, id_cond_2);
//     KRATOS_ERROR_IF(UV.size()<2) << "FitUV_BetweenSkinConditions: not enough points.\n";

//     // Ensure endpoints are exactly the first and last
//     array_1d<double,3> UV0 = UV.front();
//     array_1d<double,3> UV2 = UV.back();
//     auto t = ChordLengthParams01_UV(UV);
//     array_1d<double,3> UV1 = SolveUV_P1_LeastSquares(UV0, UV2, UV, t);

//     return MakeQuadraticNurbsUV(UV0, UV1, UV2);
// }

// /** Build a quadratic UV curve from a generic 3D edge (sample→project→LS in UV). */
// typename NurbsCurveGeometryType::Pointer MakeUV_From3DEdge(
//     const GeometryType& rEdge3D,
//     const NurbsSurfaceType& rSurface,
//     const std::size_t n_samples = 9) const
// {
//     auto UV = SampleEdgeToUV(rEdge3D, rSurface, n_samples);
//     array_1d<double,3> UV0 = UV.front();
//     array_1d<double,3> UV2 = UV.back();
//     auto t = ChordLengthParams01_UV(UV);
//     array_1d<double,3> UV1 = SolveUV_P1_LeastSquares(UV0, UV2, UV, t);
//     return MakeQuadraticNurbsUV(UV0, UV1, UV2);
// }

// /** Evaluate UV on a UV curve at param s in [0,1]. */
// inline array_1d<double,3> UV_on_curve(const NurbsCurveGeometryType& rC, const double s) const
// {
//     CoordinatesArrayType loc(3,0.0); loc[0]=s;
//     array_1d<double,3> uv; uv.clear(); uv.resize(3,false);
//     rC.GlobalCoordinates(uv, loc); // uv[0]=u, uv[1]=v
//     uv[2]=0.0;
//     return uv;
// }

// /** Coons UV mapping: UV(ξ,η) from four UV boundary curves (quadratic NURBS). */
// array_1d<double,3> CoonsUV(
//     const double xi, const double eta,
//     const NurbsCurveGeometryType& B0,  // bottom  (η=0)  : UV(ξ,0)
//     const NurbsCurveGeometryType& L0,  // left    (ξ=0)  : UV(0,η)
//     const NurbsCurveGeometryType& L1,  // right   (ξ=1)  : UV(1,η)
//     const NurbsCurveGeometryType& B1)  // top     (η=1)  : UV(ξ,1)
//     const
// {
//     const array_1d<double,3> b0 = UV_on_curve(B0, xi);
//     const array_1d<double,3> b1 = UV_on_curve(B1, xi);
//     const array_1d<double,3> l0 = UV_on_curve(L0, eta);
//     const array_1d<double,3> l1 = UV_on_curve(L1, eta);

//     const array_1d<double,3> P00 = UV_on_curve(L0, 0.0);
//     const array_1d<double,3> P01 = UV_on_curve(L0, 1.0);
//     const array_1d<double,3> P10 = UV_on_curve(L1, 0.0);
//     const array_1d<double,3> P11 = UV_on_curve(L1, 1.0);

//     const double om_x = 1.0 - xi;
//     const double om_e = 1.0 - eta;

//     // Transfinite interpolation in UV
//     array_1d<double,3> UV =
//           om_x * l0 + xi * l1
//         + om_e * b0 + eta * b1
//         - om_x * om_e * P00
//         - xi   * om_e * P10
//         - xi   * eta  * P11
//         - om_x * eta  * P01;

//     UV[2]=0.0;
//     return UV;
// }

// /** FD derivative of the UV-Coons map: ∂(u,v)/∂ξ or ∂(u,v)/∂η. */
// array_1d<double,3> CoonsUV_DerivativeFD(
//     const double xi, const double eta, const bool wrtXi,
//     const NurbsCurveGeometryType& B0,
//     const NurbsCurveGeometryType& L0,
//     const NurbsCurveGeometryType& L1,
//     const NurbsCurveGeometryType& B1,
//     const double h = 1e-6) const
// {
//     const double s1 = std::max(0.0, std::min(1.0, (wrtXi? xi+h : xi)));
//     const double s2 = std::max(0.0, std::min(1.0, (wrtXi? xi-h : xi)));
//     const double t1 = std::max(0.0, std::min(1.0, (wrtXi? eta   : eta+h)));
//     const double t2 = std::max(0.0, std::min(1.0, (wrtXi? eta   : eta-h)));

//     const array_1d<double,3> UVp = CoonsUV(s1,t1, B0,L0,L1,B1);
//     const array_1d<double,3> UVm = CoonsUV(s2,t2, B0,L0,L1,B1);
//     array_1d<double,3> d = (UVp - UVm) * (0.5 / h);
//     d[2]=0.0; return d;
// }

// /** Gauss on [0,1] from your existing helper GaussLegendreOnUnitInterval(Order,...) */

// // --- Main: build IntegrationPoint<2>(u,v, w_uv = w_ξ w_η |det J_uv|) ---
// IntegrationPointsArrayType CreateCoonsPatchGaussPointsUV(
//     const std::size_t Order,
//     const NurbsCurveGeometryType& B0,
//     const NurbsCurveGeometryType& L0,
//     const NurbsCurveGeometryType& L1,
//     const NurbsCurveGeometryType& B1) const
// {
//     IntegrationPointsArrayType gp_list;
//     gp_list.reserve(Order*Order);

//     std::vector<double> xi, w; GaussLegendreOnUnitInterval(Order, xi, w);

//     for (std::size_t i=0; i<Order; ++i)
//     for (std::size_t j=0; j<Order; ++j)
//     {
//         const double x  = xi[i];
//         const double e  = xi[j];
//         const double wi = w[i]*w[j];

//         // UV at (x,e)
//         const array_1d<double,3> UV = CoonsUV(x,e, B0,L0,L1,B1);

//         // ∂(u,v)/∂ξ and ∂(u,v)/∂η (2x2 Jacobian)
//         const array_1d<double,3> dXi  = CoonsUV_DerivativeFD(x,e,true , B0,L0,L1,B1);
//         const array_1d<double,3> dEta = CoonsUV_DerivativeFD(x,e,false, B0,L0,L1,B1);

//         // det J_uv = du/dξ * dv/dη - du/dη * dv/dξ
//         const double du_dxi  = dXi[0],  dv_dxi  = dXi[1];
//         const double du_deta = dEta[0], dv_deta = dEta[1];
//         const double detJ = std::abs(du_dxi*dv_deta - du_deta*dv_dxi);

//         // Store local UV + weight_uv
//         gp_list.emplace_back( IntegrationPoint<2>( UV[0], UV[1], wi * detJ ) );
//     }
//     return gp_list;
// }

// inline void EvalUVCurveAndTangent(
//     const NurbsCurveGeometryType& rCurve,
//     const double s,                        // s in [0,1]
//     array_1d<double,3>& rUV,              // [u,v,0]
//     array_1d<double,3>& rUV_tangent) const// [du/ds, dv/ds, 0]
// {
//     CoordinatesArrayType loc(3,0.0);
//     loc[0] = std::min(1.0, std::max(0.0, s));

//     rCurve.GlobalCoordinates(rUV, loc);
//     std::vector<array_1d<double,3>> deriv(2, ZeroVector(3));
//     rCurve.GlobalSpaceDerivatives(deriv, loc, 1); // deriv[1] = first derivative wrt s
//     rUV_tangent = deriv[1];

//     // Enforce 2D (z=0)
//     rUV[2] = 0.0;
//     rUV_tangent[2] = 0.0;
// }

// IntegrationPointsArrayType
// CreateCoonsPatchGaussPointsUV_Analytic(
//     const std::size_t Order,
//     const NurbsCurveGeometryType& rB0, // UV curve along η=0    : B0(ξ)
//     const NurbsCurveGeometryType& rL0, // UV curve along ξ=0    : L0(η)
//     const NurbsCurveGeometryType& rL1, // UV curve along ξ=1    : L1(η)
//     const NurbsCurveGeometryType& rB1) // UV curve along η=1    : B1(ξ)
//     const
// {
//     IntegrationPointsArrayType gp_list;
//     gp_list.reserve(Order*Order);

//     // Precompute Gauss nodes/weights on [0,1]
//     std::vector<double> xi, w;
//     GaussLegendreOnUnitInterval(Order, xi, w);

//     // Corner UV points (for the bilinear correction)
//     array_1d<double,3> P00, P01, P10, P11, tmp;
//     array_1d<double,3> ttmp;
//     EvalUVCurveAndTangent(rL0, 0.0, P00, tmp); // (ξ=0,η=0)
//     EvalUVCurveAndTangent(rL0, 1.0, P01, tmp); // (ξ=0,η=1)
//     EvalUVCurveAndTangent(rL1, 0.0, P10, tmp); // (ξ=1,η=0)
//     EvalUVCurveAndTangent(rL1, 1.0, P11, tmp); // (ξ=1,η=1)

//     for (std::size_t i=0; i<Order; ++i)
//     for (std::size_t j=0; j<Order; ++j)
//     {
//         const double X  = xi[i];      // ξ
//         const double E  = xi[j];      // η
//         const double wi = w[i]*w[j];

//         // Boundary values and first derivatives
//         array_1d<double,3> b0, b0_x; EvalUVCurveAndTangent(rB0, X, b0, b0_x);   // B0(ξ), dB0/dξ
//         array_1d<double,3> b1, b1_x; EvalUVCurveAndTangent(rB1, X, b1, b1_x);   // B1(ξ), dB1/dξ
//         array_1d<double,3> l0, l0_e; EvalUVCurveAndTangent(rL0, E, l0, l0_e);   // L0(η), dL0/dη
//         array_1d<double,3> l1, l1_e; EvalUVCurveAndTangent(rL1, E, l1, l1_e);   // L1(η), dL1/dη

//         const double omX = 1.0 - X;
//         const double omE = 1.0 - E;

//         // Coons UV value (not strictly needed by CreateQuadraturePointGeometries, but useful)
//         array_1d<double,3> UV =
//               omX * l0 + X * l1
//             + omE * b0 + E * b1
//             - omX*omE * P00
//             - X   *omE * P10
//             - X   *E   * P11
//             - omX*E   * P01;
//         UV[2]=0.0;

//         // Analytic partials in UV
//         array_1d<double,3> UV_xi =
//               (-l0 + l1)
//             + omE * b0_x + E * b1_x
//             + omE * P00 - omE * P10 - E * P11 + E * P01;

//         array_1d<double,3> UV_eta =
//               omX * l0_e + X * l1_e
//             + (-b0 + b1)
//             + omX * P00 + X * P10 - X * P11 - omX * P01;

//         // 2×2 Jacobian determinant in UV
//         const double du_dxi  = UV_xi[0],  dv_dxi  = UV_xi[1];
//         const double du_deta = UV_eta[0], dv_deta = UV_eta[1];
//         const double detJ_uv = std::abs(du_dxi * dv_deta - du_deta * dv_dxi);

//         // Store local (u,v) and weight = wξ wη |detJ_uv|
//         gp_list.emplace_back( IntegrationPoint<2>( UV[0], UV[1], wi * detJ_uv ) );
//     }

//     return gp_list;
// }


// // ================================================================
// // Generic p-degree Bézier/NURBS fitting in UV with fixed endpoints
// // ================================================================
// //
// // Given ordered UV samples {Q_i} between two skin conditions (inclusive),
// // fit a p-degree Bézier curve with endpoints fixed to Q_0 and Q_N:
// //     C(t) = sum_{k=0}^p B_k^p(t) P_k,  t in [0,1],  weights = 1
// // unknowns are the (p-1) interior control points P_1..P_{p-1}.
// //
// // We solve LS: minimize sum_i ||Q_i - B_0 P_0 - B_p P_p - sum_{k=1}^{p-1} B_k P_k||^2
// // Assemble normal equations for U and V independently:
// //     (A^T A) X = A^T R,  where A_{i,k} = B_k^p(t_i),  k=1..p-1
// //     R_i = Q_i - B_0 P_0 - B_p P_p
// // Add a tiny ridge on the diagonal for robustness.
// //
// // Notes:
// // - Parametrization: chord-length in UV (centripetal is also possible).
// // - All weights = 1 (polynomial). To switch to rational, extend with per-node weights.
// // - Works for p >= 2 (p=1 is trivial polyline).
// //
// // ================================================================

// /** Evaluate Bernstein basis of degree p at t in [0,1]. B has size p+1. */
// inline void BernsteinBasis(const int p, const double t, std::vector<double>& B) const
// {
//     const double s = std::min(1.0, std::max(0.0, t));
//     B.assign(p+1, 0.0);

//     // De Casteljau-like stable evaluation for Bernstein basis
//     // Start with B[0]=1 and build iteratively.
//     B[0] = 1.0;
//     const double om = 1.0 - s;
//     for (int j = 1; j <= p; ++j) {
//         double saved = 0.0;
//         for (int k = 0; k < j; ++k) {
//             const double tmp = B[k];
//             B[k] = saved + om * tmp;
//             saved = s * tmp;
//         }
//         B[j] = saved;
//     }
// }

// /** Build Bézier knot vector [0...0, 1...1] with multiplicity p+1 each side. */
// inline void BuildBezierKnots(const int p, Vector& rKnots) const
// {
//     const int nkn = 2*(p+1);
//     if ((int)rKnots.size() != nkn) rKnots.resize(nkn);
//     for (int i=0; i<=p;   ++i) rKnots[i]      = 0.0;
//     for (int i=p+1;i<nkn;++i) rKnots[i]      = 1.0;
// }

// /** Tiny dense solver for symmetric positive definite (normal eq.) with ridge. */
// inline bool SolveSPD(std::vector<std::vector<double>>& A, std::vector<double>& b) const
// {
//     // Cholesky (naive) with simple pivot guard
//     const int n = (int)A.size();
//     for (int i=0; i<n; ++i) {
//         // Diagonal
//         double sum = A[i][i];
//         for (int k=0; k<i; ++k) sum -= A[i][k]*A[i][k];
//         if (sum <= 1e-20) return false;
//         const double Lii = std::sqrt(sum);
//         A[i][i] = Lii;

//         // Off-diagonals
//         for (int j=i+1; j<n; ++j) {
//             double s = A[j][i];
//             for (int k=0; k<i; ++k) s -= A[j][k]*A[i][k];
//             A[j][i] = s / Lii;
//         }
//     }

//     // Solve L y = b
//     for (int i=0; i<n; ++i) {
//         double s = b[i];
//         for (int k=0; k<i; ++k) s -= A[i][k]*b[k];
//         b[i] = s / A[i][i];
//     }
//     // Solve L^T x = y
//     for (int i=n-1; i>=0; --i) {
//         double s = b[i];
//         for (int k=i+1; k<n; ++k) s -= A[k][i]*b[k];
//         b[i] = s / A[i][i];
//     }
//     return true;
// }

// /** Build a UV Bézier/NURBS of degree p from control points (weights=1). */
// typename NurbsCurveGeometryType::Pointer MakeBezierUV_FromControls(
//     const std::vector<array_1d<double,3>>& CtrlsUV, const int p) const
// {
//     const int ncp = p + 1;
//     KRATOS_ERROR_IF((int)CtrlsUV.size() != ncp)
//         << "MakeBezierUV_FromControls: expected " << ncp << " control points.\n";

//     PointerVector<Node> ctrl;
//     ctrl.reserve(ncp);
//     for (int i=0; i<ncp; ++i) {
//         array_1d<double,3> P = CtrlsUV[i];
//         P[2] = 0.0; // UV is 2D
//         // Virtual node (not stored in a ModelPart)
//         ctrl.push_back( Node::Pointer(new Node(1, P)) );
//     }

//     Vector knots; BuildBezierKnots(p, knots);
//     Vector w(ncp); for (int i=0;i<ncp;++i) w[i]=1.0;

//     return typename NurbsCurveGeometryType::Pointer(
//         new NurbsCurveGeometryType(ctrl, p, knots, w));
// }

// /** Generic LS fit of degree p for UV samples with fixed endpoints. */
// typename NurbsCurveGeometryType::Pointer FitBezierUV_LS_Generic(
//     const std::vector<array_1d<double,3>>& UVsamples, // ordered, inclusive
//     const int p,
//     const double ridge = 1e-12) const
// {
//     KRATOS_ERROR_IF(p < 2) << "FitBezierUV_LS_Generic: degree p must be >= 2.\n";
//     const std::size_t N = UVsamples.size();
//     KRATOS_ERROR_IF(N < (std::size_t)(p+1))
//         << "FitBezierUV_LS_Generic: need at least p+1 samples.\n";

//     // Endpoints fixed
//     const array_1d<double,3> P0 = UVsamples.front();
//     const array_1d<double,3> Pp = UVsamples.back();

//     // Parameterization
//     auto t = ChordLengthParams01_UV(UVsamples);

//     // Assemble normal equations for interior controls P1..P_{p-1}
//     const int m  = p - 1;              // number of unknown control points
//     std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0)); // A^T A
//     std::vector<double> bu(m, 0.0), bv(m, 0.0);                         // A^T R (for U and V)

//     std::vector<double> B; B.reserve(p+1);

//     for (std::size_t i=0; i<N; ++i) {
//         BernsteinBasis(p, t[i], B);

//         // Residual after removing endpoint contributions
//         const double Ru = UVsamples[i][0] - (B[0]*P0[0] + B[p]*Pp[0]);
//         const double Rv = UVsamples[i][1] - (B[0]*P0[1] + B[p]*Pp[1]);

//         // Fill (A^T A) and (A^T R)
//         for (int a=1; a<=p-1; ++a) {
//             const double Ba = B[a];
//             bu[a-1] += Ba * Ru;
//             bv[a-1] += Ba * Rv;
//             for (int b=1; b<=p-1; ++b) {
//                 M[a-1][b-1] += Ba * B[b];
//             }
//         }
//     }

//     // Ridge regularization
//     for (int d=0; d<m; ++d) M[d][d] += ridge;

//     // Solve for U and V components
//     auto Mu = M; auto Mv = M; // same matrix for both
//     bool ok_u = SolveSPD(Mu, bu);
//     bool ok_v = SolveSPD(Mv, bv);
//     KRATOS_ERROR_IF(!(ok_u && ok_v)) << "FitBezierUV_LS_Generic: SPD solve failed.\n";

//     // Collect control points: P0, P1..P_{p-1}, Pp
//     std::vector<array_1d<double,3>> Ctrls(p+1, ZeroVector(3));
//     Ctrls[0]     = P0;
//     Ctrls[p]     = Pp;
//     for (int k=1; k<=p-1; ++k) {
//         Ctrls[k][0] = bu[k-1];
//         Ctrls[k][1] = bv[k-1];
//         Ctrls[k][2] = 0.0;
//     }

//     return MakeBezierUV_FromControls(Ctrls, p);
// }

// /** Reverse orientation of a UV Bézier/NURBS (any degree): returns the same geometry, reversed. */
// typename NurbsCurveGeometryType::Pointer ReverseBezierUV_Generic(
//     const typename NurbsCurveGeometryType::Pointer& p_forward) const
// {
//     KRATOS_ERROR_IF(!p_forward) << "ReverseBezierUV_Generic: null input.\n";

//     const int p = (int)p_forward->PolynomialDegree(0);
//     const int n = (int)p_forward->size();
//     KRATOS_ERROR_IF(n != p+1) << "ReverseBezierUV_Generic: unexpected control count.\n";

//     PointerVector<Node> ctrl_rev; ctrl_rev.reserve(n);
//     for (int i=n-1; i>=0; --i) ctrl_rev.push_back(p_forward->pGetPoint(i));

//     Vector knots; BuildBezierKnots(p, knots);
//     Vector w(n); for (int i=0;i<n;++i) w[i]=1.0;

//     return typename NurbsCurveGeometryType::Pointer(
//         new NurbsCurveGeometryType(ctrl_rev, p, knots, w));
// }

// /** Fit between two skin conditions (inclusive) with general degree p in UV. */
// template <bool TIsInnerLoop>
// typename NurbsCurveGeometryType::Pointer FitUV_BetweenSkinNodes_Generic(
//     const ModelPart& rSkinSubModelPart,
//     const NurbsSurfaceType& rSurface,
//     IndexType id_node_1,
//     IndexType id_node_2,
//     const int p,
//     const double ridge = 1e-12) const
// {
//     // Reuse your robust UV collector (ordered with wrap-around).
//     auto UV = CollectSkinUVBetween<TIsInnerLoop>(rSkinSubModelPart, rSurface, id_node_1, id_node_2);
//     KRATOS_ERROR_IF(UV.size() < (std::size_t)(p+1))
//         << "FitUV_BetweenSkinConditions_Generic: not enough samples for degree p=" << p << ".\n";

//     // Ensure endpoints match exactly the first and last
//     UV.front() = UV.front();
//     UV.back()  = UV.back();

//     return FitBezierUV_LS_Generic(UV, p, ridge);
// }


// inline double Orientation(
//     const Node& p, const Node& q, const Node& r)
// {
//     return (q.X() - p.X()) * (r.Y() - p.Y()) -
//            (q.Y() - p.Y()) * (r.X() - p.X());
// }

// inline bool OnSegment(
//     const Node& p, const Node& q, const Node& r)
// {
//     return (std::min(p.X(), r.X()) <= q.X() && q.X() <= std::max(p.X(), r.X()) &&
//             std::min(p.Y(), r.Y()) <= q.Y() && q.Y() <= std::max(p.Y(), r.Y()));
// }

// inline bool SegmentsIntersect(
//     const Node& A, const Node& B, const Node& C, const Node& D)
// {
//     double o1 = Orientation(A, B, C);
//     double o2 = Orientation(A, B, D);
//     double o3 = Orientation(C, D, A);
//     double o4 = Orientation(C, D, B);

//     // Caso generale
//     if (o1 * o2 < 0.0 && o3 * o4 < 0.0)
//         return true;

//     // Collinearità (con bounding box)
//     if (std::abs(o1) < 1e-14 && OnSegment(A, C, B)) return true;
//     if (std::abs(o2) < 1e-14 && OnSegment(A, D, B)) return true;
//     if (std::abs(o3) < 1e-14 && OnSegment(C, A, D)) return true;
//     if (std::abs(o4) < 1e-14 && OnSegment(C, B, D)) return true;

//     return false;
// }

// IndexType FindClosestNodeInLayer(
//     const DynamicBinsPointerType& rStartPoint,
//     BinSearchParameters& rSearchParameters,
//     const std::string& rLayer,
//     const ModelPart& rSkinSubModelPart);

// IndexType FindClosestNodeInLayerWithDirection(
//     const DynamicBinsPointerType& rStartPoint,
//     BinSearchParameters& rSearchParameters,
//     const std::string& rLayer,
//     const ModelPart& rSkinSubModelPart,
//     const Vector& rTangentDirection);

        
// }; // Class SnakeCutSbmProcess

// }  // namespace Kratos.
