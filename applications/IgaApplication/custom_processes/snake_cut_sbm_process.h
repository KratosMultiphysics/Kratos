//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//

#pragma once

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"
#include "processes/process.h"
#include "geometries/nurbs_curve_geometry.h"
#include "snake_sbm_process.h"
#include "custom_utilities/create_breps_sbm_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class SnakeCutSbmProcess
 * @brief Process class for implementing Snake-based Surrogate Boundary Method (SBM).
 * This class provides various functions to create and manipulate surrogate boundaries
 * for Iga models in Kratos.
 */
class KRATOS_API(IGA_APPLICATION) SnakeCutSbmProcess
    : public SnakeSbmProcess
{

public:

    typedef Node NodeType;
    typedef PointerVector<Node> ContainerNodeType;
    typedef PointerVector<Point> ContainerEmbeddedNodeType;
    typedef BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;
    typedef typename ModelPart::ElementsContainerType ElementsContainerType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;
    typedef typename Properties::Pointer PropertiesPointerType;
    typedef BrepCurve<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveType;
    
    using NurbsCurveGeometryType = NurbsCurveGeometry<2, PointerVector<Node>>;

    ///@name Type Definitions
    ///@{

    /// Pointer definition of SnakeCutSbmProcess
    KRATOS_CLASS_POINTER_DEFINITION(SnakeCutSbmProcess);

    ///@name Life Cycle
    ///@{

    /// Constructor
    SnakeCutSbmProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~SnakeCutSbmProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {   
        CreateSbmExtendedGeometries();
    };
    
    void ExecuteInitialize() override
    {
        SnakeSbmProcess::CreateTheSnakeCoordinates();
    };

    void ExecuteInitializeSolutionStep() override
    {};

    void ExecuteBeforeSolutionLoop() override
    {};

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "",
            "skin_model_part_inner_initial_name" : "SkinModelPartInnerInitial",
            "skin_model_part_outer_initial_name" : "SkinModelPartOuterInitial",
            "skin_model_part_name" : "SkinModelPart",
            "echo_level" : 0,
            "lambda_inner" : 0.5,
            "lambda_outer" : 0.5,
            "number_of_inner_loops": 0,
            "cut_element_name": "",
            "cut_interface_condition_name": ""
        })" );

        return default_parameters;
    }
    
private:

    ModelPart* mpCutElementsSubModelPart = nullptr; 
    ModelPart* mpCutInterfaceSubModelPart = nullptr; 
    std::string mCutElementName;
    std::string mCutInterfaceConditionName;

    void CreateSbmExtendedGeometries();

    void FindClosestTruePointToSurrogateVertex();

    /**
     * @brief 
     * 
     * @tparam TIsInnerLoop 
     * @param rSkinSubModelPart 
     * @param rSurrogateSubModelPart 
     */
    template <bool TIsInnerLoop>
    void FindClosestTruePointToSurrogateVertex(
        const ModelPart& rSkinSubModelPart,
        const ModelPart& rSurrogateSubModelPart);

    void FindClosestTruePointToSurrogateVertexByNurbs();

    bool ProjectToSkinBoundary(
        const ModelPart* pSkinModelPart,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rProjectedPoint,
        CoordinatesArrayType& rProjectedPointLocal,
        std::vector<array_1d<double, 3>>& rCurveDerivatives,
        int nInitialGuesses);

    
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        const Vector KnotSpanSizes,
        const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const;

    /// Creates elements from geometries
    void CreateElements(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const;

    static typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer CreateBrepCurve(
        const Node::Pointer pFirstBrepPoint,
        const Node::Pointer pSecondBrepPoint, 
        const Vector& rActiveRangeKnotVector)
        {
            // Create the data for the trimming curves
            PointerVector<Node> control_points;
            control_points.push_back(pFirstBrepPoint);
            control_points.push_back(pSecondBrepPoint);
            const int polynomial_degree = 1;
            Vector knot_vector = ZeroVector(4) ;
            knot_vector[0] = rActiveRangeKnotVector[0] ;
            knot_vector[1] = rActiveRangeKnotVector[0] ;
            knot_vector[2] = rActiveRangeKnotVector[1] ;
            knot_vector[3] = rActiveRangeKnotVector[1] ;
            // Create the trimming curves
            typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_trimming_curve(
                new NurbsCurveGeometry<3, PointerVector<Node>>(
                    control_points,
                    polynomial_degree,
                    knot_vector));   
            return p_trimming_curve;
        }
    
    /// Return n×n Gauss points mapped on the curved quadrilateral
    IntegrationPointsArrayType CreateCoonsPatchGaussPoints(
        std::size_t                  Order,
        const GeometryType&         rB0,
        const GeometryType&         rL0,
        const GeometryType&         rL1,
        const GeometryType&         rB1,
        const array_1d<double,3>&    rP00,
        const array_1d<double,3>&    rP01,
        const array_1d<double,3>&    rP10,
        const array_1d<double,3>&    rP11) const;

    // --- static utility helpers ----------------------------------

    /// 1-D Gauss-Legendre nodes & weights mapped to [0,1]
    static void GaussLegendreOnUnitInterval(
        std::size_t              Order,
        std::vector<double>&     rXi,
        std::vector<double>&     rWeight);

    /// Global coordinates of a point on a Brep-curve at parameter t
    static array_1d<double,3> GlobalPoint(
        const GeometryType& rCurve,
        double               T);

    /// Coons patch mapping X(ξ,η)
    static array_1d<double,3> CoonsPoint(
        double                     Xi,
        double                     Eta,
        const GeometryType&       rB0,
        const GeometryType&       rL0,
        const GeometryType&       rL1,
        const GeometryType&       rB1,
        const array_1d<double,3>&  rP00,
        const array_1d<double,3>&  rP01,
        const array_1d<double,3>&  rP10,
        const array_1d<double,3>&  rP11);

    /// Finite-difference derivative of CoonsPoint
    static array_1d<double,3> CoonsDerivativeFD(
        double                     Xi,
        double                     Eta,
        bool                       WithRespectToXi,
        const GeometryType&       rB0,
        const GeometryType&       rL0,
        const GeometryType&       rL1,
        const GeometryType&       rB1,
        const array_1d<double,3>&  rP00,
        const array_1d<double,3>&  rP01,
        const array_1d<double,3>&  rP10,
        const array_1d<double,3>&  rP11,
        double                     Step = 1.0e-8);


    /** @brief Build control-points, knot-vector and weights of a quadratic
     *         NURBS parabola through three *existing* nodes.
     *
     *  @param pNode0  start   node  (P0)
     *  @param pNodeM  middle  node  (Pm) – parametric coord tm (default 0.5)
     *  @param pNode2  end     node  (P2)
     *  @param rCtrlPtsPointerVector  OUT pointer-vector with 3 control nodes
     *                                (new node created for the middle CP)
     *  @param rKnots               OUT Vector size 6 : [0 0 0 1 1 1]
     *  @param rWeights             OUT Vector size 3 : [1 1 1]
     *  @param tm                   parametric coordinate of Pm in (0,1)
     */
    void BuildParabolicNurbsData(
        Node::Pointer           pNode0,
        Node::Pointer           pNodeM,
        Node::Pointer           pNode2,
        PointerVector<Node>&    rCtrlPtsPointerVector,
        Vector&                 rKnots,
        Vector&                 rWeights,
        const double            tm = 0.5);
        
}; // Class SnakeCutSbmProcess

}  // namespace Kratos.
