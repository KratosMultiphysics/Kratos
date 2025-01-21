//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_SNAKE_SBM_UTILITIES_NEW_H_INCLUDED )
#define  KRATOS_SNAKE_SBM_UTILITIES_NEW_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    class KRATOS_API(IGA_APPLICATION) SnakeSBMUtilities
    {
    public:
        typedef std::size_t IndexType;
        typedef std::size_t SizeType;
        typedef Node NodeType;

        typedef Geometry<NodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;
        using NurbsCurveGeometryPointerType = NurbsCurveGeometry<2, PointerVector<Node>>::Pointer;
        using CoordinatesArrayType = Geometry<NodeType>::CoordinatesArrayType;

        static void CreateTheSnakeCoordinates(ModelPart& iga_model_part, 
                                            ModelPart& skin_model_part_inner_initial,
                                            ModelPart& skin_model_part_outer_initial,
                                            ModelPart& skin_model_part, 
                                            int rEchoLevel, 
                                            Vector& knot_vector_u, 
                                            Vector& knot_vector_v, 
                                            const Parameters mParameters);
        
    
        static void CreateTheSnakeCoordinates(ModelPart& iga_model_part, 
                                            ModelPart& skin_model_part_initial,
                                            ModelPart& skin_model_part, 
                                            int rEchoLevel, 
                                            Vector& knot_vector_u, 
                                            Vector& knot_vector_v, 
                                            const Parameters mParameters,
                                            bool is_inner_loop); // FIXME : 

                                              
        static void SnakeStep(ModelPart& skin_model_part, 
                            std::vector<std::vector<std::vector<int>>> &knot_spans_available, 
                            int idMatrixKnotSpansAvailable, 
                            std::vector<std::vector<int>> knot_spans_uv, 
                            std::vector<std::vector<double>> xy_coord_i_cond, 
                            Vector knot_step_uv, 
                            Vector starting_pos_uv);
        
        static void SnakeStepNurbs(ModelPart& skin_model_part, 
                            std::vector<std::vector<std::vector<int>>> &knot_spans_available, 
                            int idMatrix, 
                            std::vector<std::vector<int>> knot_spans_uv, 
                            std::vector<std::vector<double>> xy_coord, 
                            Vector knot_step_uv, 
                            const Vector starting_pos_uv,
                            const std::vector<double> local_coords,
                            const NurbsCurveGeometryPointerType p_curve);
        
        static void SnakeStep3D(ModelPart& skin_model_part, std::vector<std::vector<std::vector<std::vector<int>>>> &knot_spans_available, int idMatrixKnotSpansAvailable, 
                            int knot_span_u_1st_point, int knot_span_u_2nd_point, int knot_span_u_3rd_point, 
                            int knot_span_v_1st_point, int knot_span_v_2nd_point, int knot_span_v_3rd_point, 
                            int knot_span_w_1st_point, int knot_span_w_2nd_point, int knot_span_w_3rd_point, 
                            double& x_true_boundary1, double& x_true_boundary2, double& x_true_boundary3, 
                            double& y_true_boundary1, double& y_true_boundary2, double& y_true_boundary3, 
                            double& z_true_boundary1, double& z_true_boundary2, double& z_true_boundary3, 
                            double& knot_step_u, double& knot_step_v, double& knot_step_z,  array_1d<IndexType, 3>& ordered_ids);
    private:
        ///@name Iga functionalities
        ///@{
        using PointType = Node;
        using PointTypePointer = Node::Pointer;
        using PointVector = std::vector<PointType::Pointer>;
        using PointIterator = std::vector<PointType::Pointer>::iterator;
        using DistanceVector = std::vector<double>;
        using DistanceIterator = std::vector<double>::iterator;
        using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
        using PointerType = DynamicBins::PointerType;


        static bool isPointInsideSkinBoundary(Point& pointToSearch, DynamicBins& testBins, ModelPart& skin_model_part_in);

        static void MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>> & knot_spans_available, 
                                            int idMatrix,DynamicBins& testBin, 
                                            ModelPart& skin_model_part, 
                                            double lambda, 
                                            std::vector<int>& n_knot_spans_uv, 
                                            Vector& knot_step_uv, 
                                            const Vector& starting_pos_uv);

        static void CreateSurrogateBuondaryFromSnake_inner (std::vector<std::vector<std::vector<int>>> & knot_spans_available, 
                                                            int idMatrix, 
                                                            ModelPart& surrogate_model_part_inner, 
                                                            std::vector<int>& n_knot_spans_uv, 
                                                            Vector& knot_vector_u,
                                                            Vector& knot_vector_v,
                                                            const Vector& starting_pos_uv);

        static void CreateSurrogateBuondaryFromSnake_outer (DynamicBins& testBin_out, 
                                                            ModelPart& initial_skin_model_part_out, 
                                                            std::vector<std::vector<std::vector<int>>> & knot_spans_available, 
                                                            int idMatrix, 
                                                            ModelPart& surrogate_model_part_outer, 
                                                            std::vector<int>& n_knot_spans_uv, 
                                                            Vector& knot_vector_u, 
                                                            Vector&  knot_vector_v,
                                                            const Vector& starting_pos_uv);

    }; // Class DirectorUtilities

}  // namespace Kratos.

#endif // KRATOS_DIRECTOR_UTILITIES_NEW_H_INCLUDED  defined
