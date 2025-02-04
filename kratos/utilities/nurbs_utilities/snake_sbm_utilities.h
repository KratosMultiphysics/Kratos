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
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    /**
     * @class SnakeSbmUtilities
     * @brief Utility class for implementing Snake-based Surrogate Boundary Method (SBM).
     * This class provides various functions to create and manipulate snake-like surrogate boundaries
     * for finite element models in Kratos.
     */
    class KRATOS_API(IGA_APPLICATION) SnakeSbmUtilities
    {
    public:
        /**
         * @brief Creates the initial snake coordinates for 2D or 3D skin.
         * @param rIgaModelPart Model part containing the IGA geometry.
         * @param rSkinModelPartInnerInitial Model part for the initial inner skin boundary.
         * @param rSkinModelPartOuterInitial Model part for the initial outer skin boundary.
         * @param rSkinModelPart Model part where the snake coordinates will be stored.
         * @param rEchoLevel Verbosity level for logging.
         * @param knotVectorU Knot vector in the U-direction.
         * @param knotVectorV Knot vector in the V-direction.
         * @param mParameters Input parameters for the process 
         *                      ex:
         *                      {
                                    "sbm_parameters": {
                                        "lambda_inner": 0.5,
                                        "lambda_outer": 0.5,
                                        "number_of_inner_loops": 1
                                    }
                                }
         */
        static void CreateTheSnakeCoordinates(ModelPart& rIgaModelPart, 
                                            ModelPart& rSkinModelPartInnerInitial,
                                            ModelPart& rSkinModelPartOuterInitial,
                                            ModelPart& rSkinModelPart, 
                                            int rEchoLevel, 
                                            Vector& knotVectorU, 
                                            Vector& knotVectorV, 
                                            const Parameters mParameters);
        
        /**
         * @brief Creates the snake coordinates for a specific loop (inner or outer).
         * @param rIgaModelPart Model part containing the IGA geometry.
         * @param rSkinModelPartInitial Model part for the initial skin boundary.
         * @param rSkinModelPart Model part where the snake coordinates will be stored.
         * @param rEchoLevel Verbosity level for logging.
         * @param knotVectorU Knot vector in the U-direction.
         * @param knotVectorV Knot vector in the V-direction.
         * @param mParameters Input parameters for the process.
         *                      ex:
         *                      {
                                    "sbm_parameters": {
                                        "lambda_inner": 0.5,
                                        "lambda_outer": 0.5,
                                        "number_of_inner_loops": 1
                                    }
                                }
         * @param isInnerLoop Boolean flag to indicate if the loop is inner or outer.
         */
        static void CreateTheSnakeCoordinates(ModelPart& rIgaModelPart, 
                                            ModelPart& rSkinModelPartInitial,
                                            ModelPart& rSkinModelPart, 
                                            int rEchoLevel, 
                                            Vector& knotVectorU, 
                                            Vector& knotVectorV, 
                                            const Parameters mParameters,
                                            bool isInnerLoop);

        /**
         * @brief Performs a single step in the snake algorithm for 2D models.
         * @param rSkinModelPart The skin model part to be updated.
         * @param knotSpansAvailable Knot spans available for the snake.
         * @param idMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
         * @param knotSpansUV Knot spans in UV coordinates.
         * @param conditionCoords XY coordinates of control points.
         * @param knotStepUV Step size in UV space.
         * @param startingPositionUV Starting position in UV space.
         */
        static void SnakeStep(ModelPart& rSkinModelPart, 
                            std::vector<std::vector<std::vector<int>>> &knotSpansAvailable, 
                            int idMatrixKnotSpansAvailable, 
                            std::vector<std::vector<int>> knotSpansUV, 
                            std::vector<std::vector<double>> conditionCoords, 
                            Vector knotStepUV, 
                            Vector startingPositionUV);

    private:
        ///@name IGA functionalities
        ///@{
        using PointType = Node;
        using PointTypePointer = Node::Pointer;
        using PointVector = std::vector<PointType::Pointer>;
        using PointIterator = std::vector<PointType::Pointer>::iterator;
        using DistanceVector = std::vector<double>;
        using DistanceIterator = std::vector<double>::iterator;
        using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using PointerType = DynamicBins::PointerType;

        /**
         * @brief Checks if a point is inside the given skin boundary.
         * @param pointToSearch The point to be checked.
         * @param testBins Spatial bins to search for nearby points.
         * @param rSkinModelPartIn Input skin model part.
         * @return True if the point is inside the boundary, false otherwise.
         */
        static bool IsPointInsideSkinBoundary(Point& pointToSearch, 
                                              DynamicBins& testBins,
                                              ModelPart& rSkinModelPartIn);

        /**
         * @brief Marks the knot spans as available for the snake algorithm.
         * @param knotSpansAvailable Reference to the knot spans availability matrix.
         * @param idMatrix ID of the matrix being updated.
         * @param testBin Spatial bins to search for nearby points.
         * @param rSkinModelPart Skin model part being analyzed.
         * @param lambda Relaxation parameter.
         * @param nKnotSpansUV Number of knot spans in UV directions.
         * @param knotStepUV Step size in UV space.
         * @param startingPositionUV Starting position in UV space.
         */
        static void MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, 
                                            int idMatrix,
                                            DynamicBins& testBin, 
                                            ModelPart& rSkinModelPart, 
                                            double lambda, 
                                            std::vector<int>& nKnotSpansUV, 
                                            array_1d<double, 2>& knotStepUV, 
                                            const Vector& startingPositionUV);

        /**
         * @brief Creates the inner surrogate boundary from the snake algorithm.
         * @param knotSpansAvailable Reference to the knot spans availability matrix.
         * @param idMatrix ID of the matrix being updated.
         * @param rSurrogateModelPartInner Model part representing the inner boundary.
         * @param nKnotSpansUV Number of knot spans in UV directions.
         * @param knotVectorU Knot vector in the U-direction.
         * @param knotVectorV Knot vector in the V-direction.
         * @param startingPositionUV Starting position in UV space.
         */
        static void CreateSurrogateBuondaryFromSnakeInner (std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, 
                                                            int idMatrix, 
                                                            ModelPart& rSurrogateModelPartInner, 
                                                            std::vector<int>& nKnotSpansUV, 
                                                            Vector& knotVectorU,
                                                            Vector& knotVectorV,
                                                            const Vector& startingPositionUV);

        /**
         * @brief Creates the outer surrogate boundary from the snake algorithm.
         * @param testBinOuter Spatial bins for testing.
         * @param rInitialSkinModelPartOuter Initial outer skin model part.
         * @param knotSpansAvailable Reference to the knot spans availability matrix.
         * @param idMatrix ID of the matrix being updated.
         * @param rSurrogateModelPartOuter Model part representing the outer boundary.
         * @param nKnotSpansUV Number of knot spans in UV directions.
         * @param knotVectorU Knot vector in the U-direction.
         * @param knotVectorV Knot vector in the V-direction.
         * @param startingPositionUV Starting position in UV space.
         */
        static void CreateSurrogateBuondaryFromSnakeOuter (DynamicBins& testBinOuter, 
                                                            ModelPart& rInitialSkinModelPartOuter, 
                                                            std::vector<std::vector<std::vector<int>>> & knotSpansAvailable, 
                                                            int idMatrix, 
                                                            ModelPart& rSurrogateModelPartOuter, 
                                                            std::vector<int>& nKnotSpansUV, 
                                                            Vector& knotVectorU, 
                                                            Vector&  knotVectorV,
                                                            const Vector& startingPositionUV);

    }; // Class SnakeSbmUtilities

}  // namespace Kratos.
