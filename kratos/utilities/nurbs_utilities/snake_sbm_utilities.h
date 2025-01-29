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
     * @class SnakeSBMUtilities
     * @brief Utility class for implementing Snake-based Surrogate Boundary Method (SBM).
     * This class provides various functions to create and manipulate snake-like surrogate boundaries
     * for finite element models in Kratos.
     */
    class KRATOS_API(IGA_APPLICATION) SnakeSBMUtilities
    {
    public:
        /**
         * @brief Creates the initial snake coordinates for 2D or 3D skin.
         * @param iga_model_part Model part containing the IGA geometry.
         * @param skin_model_part_inner_initial Model part for the initial inner skin boundary.
         * @param skin_model_part_outer_initial Model part for the initial outer skin boundary.
         * @param skin_model_part Model part where the snake coordinates will be stored.
         * @param rEchoLevel Verbosity level for logging.
         * @param knot_vector_u Knot vector in the U-direction.
         * @param knot_vector_v Knot vector in the V-direction.
         * @param mParameters Input parameters for the process.
         */
        static void CreateTheSnakeCoordinates(ModelPart& iga_model_part, 
                                            ModelPart& skin_model_part_inner_initial,
                                            ModelPart& skin_model_part_outer_initial,
                                            ModelPart& skin_model_part, 
                                            int rEchoLevel, 
                                            Vector& knot_vector_u, 
                                            Vector& knot_vector_v, 
                                            const Parameters mParameters);
        
        /**
         * @brief Creates the snake coordinates for a specific loop (inner or outer).
         * @param iga_model_part Model part containing the IGA geometry.
         * @param skin_model_part_initial Model part for the initial skin boundary.
         * @param skin_model_part Model part where the snake coordinates will be stored.
         * @param rEchoLevel Verbosity level for logging.
         * @param knot_vector_u Knot vector in the U-direction.
         * @param knot_vector_v Knot vector in the V-direction.
         * @param mParameters Input parameters for the process.
         * @param is_inner_loop Boolean flag to indicate if the loop is inner or outer.
         */
        static void CreateTheSnakeCoordinates(ModelPart& iga_model_part, 
                                            ModelPart& skin_model_part_initial,
                                            ModelPart& skin_model_part, 
                                            int rEchoLevel, 
                                            Vector& knot_vector_u, 
                                            Vector& knot_vector_v, 
                                            const Parameters mParameters,
                                            bool is_inner_loop);

        /**
         * @brief Performs a single step in the snake algorithm for 2D models.
         * @param skin_model_part The skin model part to be updated.
         * @param knot_spans_available Knot spans available for the snake.
         * @param idMatrixKnotSpansAvailable ID of the matrix tracking available knot spans.
         * @param knot_spans_uv Knot spans in UV coordinates.
         * @param xy_coord_i_cond XY coordinates of control points.
         * @param knot_step_uv Step size in UV space.
         * @param starting_pos_uv Starting position in UV space.
         */
        static void SnakeStep(ModelPart& skin_model_part, 
                            std::vector<std::vector<std::vector<int>>> &knot_spans_available, 
                            int idMatrixKnotSpansAvailable, 
                            std::vector<std::vector<int>> knot_spans_uv, 
                            std::vector<std::vector<double>> xy_coord_i_cond, 
                            Vector knot_step_uv, 
                            Vector starting_pos_uv);

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
         * @param skin_model_part_in Input skin model part.
         * @return True if the point is inside the boundary, false otherwise.
         */
        static bool IsPointInsideSkinBoundary(Point& pointToSearch, 
                                                DynamicBins& testBins,
                                                ModelPart& skin_model_part_in);

        /**
         * @brief Marks the knot spans as available for the snake algorithm.
         * @param knot_spans_available Reference to the knot spans availability matrix.
         * @param idMatrix ID of the matrix being updated.
         * @param testBin Spatial bins to search for nearby points.
         * @param skin_model_part Skin model part being analyzed.
         * @param lambda Relaxation parameter.
         * @param n_knot_spans_uv Number of knot spans in UV directions.
         * @param knot_step_uv Step size in UV space.
         * @param starting_pos_uv Starting position in UV space.
         */
        static void MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>> & knot_spans_available, 
                                            int idMatrix,DynamicBins& testBin, 
                                            ModelPart& skin_model_part, 
                                            double lambda, 
                                            std::vector<int>& n_knot_spans_uv, 
                                            array_1d<double, 2>& knot_step_uv, 
                                            const Vector& starting_pos_uv);

        /**
         * @brief Creates the inner surrogate boundary from the snake algorithm.
         * @param knot_spans_available Reference to the knot spans availability matrix.
         * @param idMatrix ID of the matrix being updated.
         * @param surrogate_model_part_inner Model part representing the inner boundary.
         * @param n_knot_spans_uv Number of knot spans in UV directions.
         * @param knot_vector_u Knot vector in the U-direction.
         * @param knot_vector_v Knot vector in the V-direction.
         * @param starting_pos_uv Starting position in UV space.
         */
        static void CreateSurrogateBuondaryFromSnakeInner (std::vector<std::vector<std::vector<int>>> & knot_spans_available, 
                                                            int idMatrix, 
                                                            ModelPart& surrogate_model_part_inner, 
                                                            std::vector<int>& n_knot_spans_uv, 
                                                            Vector& knot_vector_u,
                                                            Vector& knot_vector_v,
                                                            const Vector& starting_pos_uv);

        /**
         * @brief Creates the outer surrogate boundary from the snake algorithm.
         * @param testBin_out Spatial bins for testing.
         * @param initial_skin_model_part_out Initial outer skin model part.
         * @param knot_spans_available Reference to the knot spans availability matrix.
         * @param idMatrix ID of the matrix being updated.
         * @param surrogate_model_part_outer Model part representing the outer boundary.
         * @param n_knot_spans_uv Number of knot spans in UV directions.
         * @param knot_vector_u Knot vector in the U-direction.
         * @param knot_vector_v Knot vector in the V-direction.
         * @param starting_pos_uv Starting position in UV space.
         */
        static void CreateSurrogateBuondaryFromSnakeOuter (DynamicBins& testBin_out, 
                                                            ModelPart& initial_skin_model_part_out, 
                                                            std::vector<std::vector<std::vector<int>>> & knot_spans_available, 
                                                            int idMatrix, 
                                                            ModelPart& surrogate_model_part_outer, 
                                                            std::vector<int>& n_knot_spans_uv, 
                                                            Vector& knot_vector_u, 
                                                            Vector&  knot_vector_v,
                                                            const Vector& starting_pos_uv);

    }; // Class SnakeSBMUtilities

}  // namespace Kratos.
