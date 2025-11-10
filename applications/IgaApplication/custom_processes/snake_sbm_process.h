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
            "number_of_inner_loops": 0
        })" );

        return default_parameters;
    }
    
private:

    Model* mpModel = nullptr;
    Parameters mThisParameters;
    IndexType mEchoLevel;
    double mLambdaInner;
    double mLambdaOuter;
    std::size_t mNumberOfInnerLoops;
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

    /**
    * @brief Creates the initial snake coordinates for 2D skin.
    */
    void CreateTheSnakeCoordinates();

    /**
     * @brief Create a The Snake Coordinates 2 D object
     * 
     */
   void CreateTheSnakeCoordinates2D();

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
        const std::size_t NumberOfInnerLoops,
        const double Lambda,
        IndexType EchoLevel,
        ModelPart& rIgaModelPart,
        ModelPart& rSkinModelPart 
        );

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
        const std::vector<std::vector<double>>& ConditionCoords, 
        const Vector KnotStepUV, 
        const Vector StartingPositionUV,
        ModelPart& rSkinModelPart, 
        std::vector<std::vector<std::vector<int>>> & rKnotSpansAvailable
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
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartInner 
     */
    static void CreateSurrogateBuondaryFromSnakeInner(
        const int IdMatrix,
        const ModelPart& SkinModelPartInner,
        DynamicBins& rPointsBinInner,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        std::vector<std::vector<std::vector<int>>> & rKnotSpansAvailable,
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
        const ModelPart& SkinModelPartOuter,
        DynamicBins& rPointsBinOuter,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rStartingPositionUV,
        std::vector<std::vector<std::vector<int>>> & rKnotSpansAvailable,
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
     * @brief Create a The Snake Coordinates 3 D object
     * 
     */
   void CreateTheSnakeCoordinates3D();


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
    static void CreateTheSnakeCoordinates3D(
        const ModelPart& rSkinModelPartInitial,
        const std::size_t NumberOfInnerLoops,
        const double Lambda,
        IndexType EchoLevel,
        ModelPart& rIgaModelPart,
        ModelPart& rSkinModelPart 
        );
    
    /**
     * @brief Checks if the knot span is at the border of the parameter sapce
     * 
     * @param rKnotSpanUV 
     * @param rNumberKnotSpansUV 
     * @return true 
     * @return false 
     */
    static bool IsInside3D(
        const std::vector<std::vector<int>> & rKnotSpanUV,
        const std::vector<int> &rNumberKnotSpansUV
        ); 

    
    /**
     * @brief Performs a single step in the snake algorithm for 3D models.
     * 
     * @param IdMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
     * @param rKnotSpansUV Knot spans in UV coordinates.
     * @param ConditionCoords XY coordinates of control points.
     * @param KnotStepUV Step size in UV space.
     * @param StartingPositionUV Starting position in UV space.
     * @param rSkinModelPart The skin model part to be updated.
     * @param rKnotSpansAvailable Knot spans available for the snake.
     */
    static void SnakeStep3D(
        const int IdMatrixKnotSpansAvailable, 
        const std::vector<std::vector<int>>& rKnotSpansUV, 
        const std::vector<std::vector<double>>& ConditionCoords, 
        const Vector KnotStepUV, 
        const Vector StartingPositionUV,
        ModelPart& rSkinModelPart, 
        std::vector<std::vector<std::vector<std::vector<int>>>> & rKnotSpansAvailable,
        array_1d<IndexType, 3>& ordered_ids
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
    static void MarkKnotSpansAvailable3D(
        const int IdMatrix,
        DynamicBins& rPointsBin, 
        const ModelPart& rSkinModelPart, 
        const double Lambda, 
        const std::vector<int>& rNumberKnotSpansUVW, 
        const array_1d<double, 3>& rKnotStepUVW, 
        const Vector& rStartingPositionUV,
        std::vector<std::vector<std::vector<std::vector<int>>>> & rKnotSpansAvailable
        );

    /**
     * @brief Checks if a point is inside the given skin boundary.
     * @param rPointToSearch The point to be checked.
     * @param rPointsBin Spatial bins to search for nearby points.
     * @param rSkinModelPartIn Input skin model part.
     * @return True if the point is inside the boundary, false otherwise.
     */
    static bool IsPointInsideSkinBoundary3D(
        const Point& rPointToSearch, 
        DynamicBins& rPointsBin,
        const ModelPart& rSkinModelPartIn
        );

    
    /**
     * @brief Create a Surrogate Buondary From Snake Inner 3 D object
     * 
     * @param IdMatrix 
     * @param SkinModelPartInner 
     * @param rPointsBinInner 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rKnotVectorW 
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartInner 
     */
    static void CreateSurrogateBuondaryFromSnakeInner3D(
        const int IdMatrix,
        const ModelPart& SkinModelPartInner,
        DynamicBins& rPointsBinInner,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rKnotVectorW,
        std::vector<std::vector<std::vector<std::vector<int>>>> & rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartInner
        );

    /**
     * @brief Create a Surrogate Buondary From Snake Outer 3 D object
     * 
     * @param IdMatrix 
     * @param SkinModelPartOuter 
     * @param rPointsBinOuter 
     * @param rNumberKnotSpans 
     * @param rKnotVectorU 
     * @param rKnotVectorV 
     * @param rKnotVectorW 
     * @param rStartingPositionUVW 
     * @param rKnotSpansAvailable 
     * @param rSurrogateModelPartOuter 
     */
    static void CreateSurrogateBuondaryFromSnakeOuter3D(
        const int IdMatrix,
        const ModelPart& SkinModelPartOuter,
        DynamicBins& rPointsBinOuter,
        const std::vector<int>& rNumberKnotSpans,
        const Vector& rKnotVectorU,
        const Vector& rKnotVectorV,
        const Vector& rKnotVectorW,
        const Vector& rStartingPositionUVW,
        std::vector<std::vector<std::vector<std::vector<int>>>> & rKnotSpansAvailable,
        ModelPart& rSurrogateModelPartOuter
        );

}; // Class SnakeSbmProcess

}  // namespace Kratos.
