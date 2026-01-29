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

#include <algorithm>
#include <cmath>

#include <queue>
#include <utility>


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class SnakeSbmProcess
 * @brief Process class for implementing Snake-based Surrogate Boundary Method (SBM).
 * This class provides various functions to create and manipulate surrogate boundaries
 * for Iga models in Kratos.
 */
class KRATOS_API(IGA_APPLICATION) SnakeSbmProcess
    : public Process
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SnakeSbmProcess
    KRATOS_CLASS_POINTER_DEFINITION(SnakeSbmProcess);

    ///@name Life Cycle
    ///@{

    /// Constructor
    SnakeSbmProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~SnakeSbmProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        CreateTheSnakeCoordinates();
    };
    
    void ExecuteInitialize() override
    {};

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
            "number_initial_points_if_importing_nurbs": 5000
        })" );

        return default_parameters;
    }

protected:
    /**
    * @brief Creates the initial snake coordinates for 2D skin.
    */
    void CreateTheSnakeCoordinates(bool RemoveIslands = false);

    /**
     * @brief Create a The Snake Coordinates object
     * 
     * @tparam TIsInnerLoop 
     * @param rIgaModelPart 
     * @param rSkinModelPartInitial 
     * @param rSkinModelPart 
     * @param NumberOfInnerLoops 
     * @param Lambda 
     */
    template <bool TIsInnerLoop>
    static void CreateTheSnakeCoordinates(
        const ModelPart& rSkinModelPartInitial,
        const std::size_t NumberOfLoops,
        const double Lambda,
        IndexType EchoLevel,
        ModelPart& rIgaModelPart,
        ModelPart& rSkinModelPart,
        const int NumberInitialPointsIfImportingNurbs,
        bool RemoveIslands = false
        );


    Model* mpModel = nullptr;
    Parameters mThisParameters;
    IndexType mEchoLevel;
    double mLambdaInner;
    double mLambdaOuter;
    std::size_t mNumberOfInnerLoops;
    int mNumberInitialPointsIfImportingNurbs;
    ModelPart* mpIgaModelPart = nullptr; 
    ModelPart* mpSkinModelPartInnerInitial = nullptr; 
    ModelPart* mpSkinModelPartOuterInitial = nullptr; 
    ModelPart* mpSkinModelPart = nullptr; 

    using PointType = Node;
    using PointTypePointer = Node::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;
    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    using DynamicBinsPointerType = DynamicBins::PointerType;
    using NurbsCurveGeometryPointerType = NurbsCurveGeometry<2, PointerVector<Node>>::Pointer;
    using CoordinatesArrayType = Geometry<PointType>::CoordinatesArrayType;

    /**
     * @brief 
     * 
     * @param rPointP 
     * @param rPointQ 
     * @param rPointR 
     * @return double 
     */
    static double Orientation(
        const Node& rPointP,
        const Node& rPointQ,
        const Node& rPointR)
    {
        return (rPointQ.X() - rPointP.X()) * (rPointR.Y() - rPointP.Y()) -
               (rPointQ.Y() - rPointP.Y()) * (rPointR.X() - rPointP.X());
    }

    /**
     * @brief 
     * 
     * @param rPointP 
     * @param rPointQ 
     * @param rPointR 
     * @return true 
     * @return false 
     */
    static bool OnSegment(
        const Node& rPointP,
        const Node& rPointQ,
        const Node& rPointR)
    {
        return (std::min(rPointP.X(), rPointR.X()) <= rPointQ.X() && rPointQ.X() <= std::max(rPointP.X(), rPointR.X()) &&
                std::min(rPointP.Y(), rPointR.Y()) <= rPointQ.Y() && rPointQ.Y() <= std::max(rPointP.Y(), rPointR.Y()));
    }

    /**
     * @brief 
     * 
     * @param rPointA 
     * @param rPointB 
     * @param rPointC 
     * @param rPointD 
     * @return true 
     * @return false 
     */
    static bool SegmentsIntersect(
        const Node& rPointA,
        const Node& rPointB,
        const Node& rPointC,
        const Node& rPointD)
    {
        const double orientation_1 = Orientation(rPointA, rPointB, rPointC);
        const double orientation_2 = Orientation(rPointA, rPointB, rPointD);
        const double orientation_3 = Orientation(rPointC, rPointD, rPointA);
        const double orientation_4 = Orientation(rPointC, rPointD, rPointB);

        if (orientation_1 * orientation_2 < 0.0 && orientation_3 * orientation_4 < 0.0) {
            return true;
        }

        if (std::abs(orientation_1) < 1e-14 && OnSegment(rPointA, rPointC, rPointB)) return true;
        if (std::abs(orientation_2) < 1e-14 && OnSegment(rPointA, rPointD, rPointB)) return true;
        if (std::abs(orientation_3) < 1e-14 && OnSegment(rPointC, rPointA, rPointD)) return true;
        if (std::abs(orientation_4) < 1e-14 && OnSegment(rPointC, rPointB, rPointD)) return true;

        return false;
    }

private:

    /**
     * @brief Performs a single step in the snake algorithm for 2D models.
     * 
     * @param IdMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
     * @param rKnotSpansUV Knot spans in UV coordinates.
     * @param ConditionCoords XY coordinates of control points.
     * @param KnotStepUV Step size in UV space.
     * @param StartingPositionUV Starting position in UV space.
     * @param rSkinModelPart The skin model part to be updated.
     * @param rKnotSpansAvailable Knot spans available for the snake.
     */
    static void SnakeStep(
        const int IdMatrixKnotSpansAvailable, 
        const std::vector<std::vector<int>>& rKnotSpansUV, 
        const std::vector<std::vector<double>>& rConditionCoordinates, 
        const array_1d<double, 2>& rKnotStepUV, 
        const array_1d<double, 2>& rStartingPositionUV,
        ModelPart& rSkinModelPart, 
        std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable
        );
    
    /**
     * @brief Performs a single step in the snake algorithm for 2D models with NURBS geometries.
     * 
     * @param IdMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
     * @param rKnotSpansUV  Knot spans in UV coordinates.
     * @param rConditionCoord  XY coordinates of control points.
     * @param rKnotStepUV  Step size in UV space.
     * @param rStartingPosition Starting position in UV space.
     * @param rLocalCoords Local Coords of the point in the NURBS curve geometry.
     * @param rpCurve Pointer of the NURBS curve geometry.
     * @param rSkinModelPart The skin model part to be updated.
     * @param rKnotSpansAvailable  Knot spans available for the snake.
     */
    static void SnakeStepNurbs(
            const int IdMatrixKnotSpansAvailable, 
            const std::vector<std::vector<int>>& rKnotSpansUV, 
            const std::vector<std::vector<double>>& rConditionCoord, 
            const array_1d<double, 2>& rKnotStepUV, 
            const array_1d<double, 2>& rStartingPosition,
            const std::vector<double>& rLocalCoords,
            const NurbsCurveGeometryPointerType& rpCurve,
            ModelPart& rSkinModelPart, 
            std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable
            );

    /**
     * @brief Checks if a point is inside the given skin boundary.
     * @param rPointToSearch The point to be checked.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPartIn Input skin model part.
     * @return True if the point is inside the boundary, false otherwise.
     */
    static bool IsPointInsideSkinBoundary(
        const Point& rPointToSearch, 
        DynamicBins& rPointsBin,
        const ModelPart& rSkinModelPartIn
        );

    /**
     * @brief Marks the knot spans as available for the snake algorithm.
     * @param IdMatrix ID of the matrix being updated.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPart Skin model part being analyzed.
     * @param Lambda Relaxation parameter.
     * @param rNumberKnotSpansUV Number of knot spans in UV directions.
     * @param rKnotStepUV Step size in UV space.
     * @param rStartingPositionUV Starting position in UV space.
     * @param rKnotSpansAvailable Reference to the knot spans availability matrix.
     */
    static void MarkKnotSpansAvailable(
        const int IdMatrix,
        DynamicBins& rPointsBin, 
        const ModelPart& rSkinModelPart, 
        const double Lambda, 
        const std::vector<int>& rNumberKnotSpansUV, 
        const array_1d<double, 2>& rKnotStepUV, 
        const Vector& rStartingPositionUV,
        std::vector<std::vector<std::vector<int>>> & rKnotSpansAvailable
        );

    /**
     * @brief Create a Surrogate Buondary From Snake Inner object
     * 
     * @param IdMatrix 
     * @param SkinModelPartInner 
     * @param rPointsBinInner 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rStartingPositionUV
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartInner 
     */
    static void CreateSurrogateBuondaryFromSnakeInner(
        const int IdMatrix,
        const ModelPart& rSkinModelPartInner,
        DynamicBins& rPointsBinInner,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rStartingPositionUV,
        std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartInner
        );

    /**
     * @brief Create a Surrogate Buondary From Snake Outer object
     * 
     * @param IdMatrix 
     * @param SkinModelPartOuter 
     * @param rPointsBinOuter 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rStartingPositionUV 
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartOuter 
     */
    static void CreateSurrogateBuondaryFromSnakeOuter(
        const int IdMatrix,
        const ModelPart& rSkinModelPartOuter,
        DynamicBins& rPointsBinOuter,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rStartingPositionUV,
        std::vector<std::vector<std::vector<int>>>& rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartOuter
        );
    
    /**
     * @brief Checks if the knot span is at the border of the parameter sapce
     * 
     * @param rKnotSpanUV 
     * @param rNumberKnotSpansUV 
     * @return true 
     * @return false 
     */
    static bool IsInside(
        const std::vector<std::vector<int>> & rKnotSpanUV,
        const std::vector<int> &rNumberKnotSpansUV
        ); 

    
    /**
     * @brief Remove the islands surrogate domain disconnected from the main one
     * 
     * @tparam TIsInnerLoop 
     * @param grid 
     */
    template <bool TIsInnerLoop>
    static void KeepLargestZeroIsland(std::vector<std::vector<int>>& rGrid);

}; // Class SnakeSbmProcess

}  // namespace Kratos.
