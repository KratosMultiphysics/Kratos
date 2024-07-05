//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_SNAKE_SBM_UTILITIES_H_INCLUDED )
#define  KRATOS_SNAKE_SBM_UTILITIES_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "spatial_containers/bins_dynamic.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    class KRATOS_API(IGA_APPLICATION) SnakeSBMUtilities
    {
    public:

        static void CreateTheSnakeCoordinates(ModelPart& iga_model_part, ModelPart& skin_model_part_in, ModelPart& skin_model_part_out, 
                                              ModelPart& initial_skin_model_part_in, ModelPart& initial_skin_model_part_out,  int rEchoLevel, 
                                              Vector& knot_vector_u_complete, Vector& knot_vector_v_complete, double& knot_step_u, double& knot_step_v, 
                                              const Parameters refinements_parameters, const Parameters mParameters, ModelPart& surrogate_model_part_inner, 
                                              ModelPart& surrogate_model_part_outer);

        static void CreateTheSnakeCoordinates3D(ModelPart& iga_model_part, ModelPart& skin_model_part_in, ModelPart& skin_model_part_out, 
                                              ModelPart& initial_skin_model_part_in, ModelPart& initial_skin_model_part_out,  int rEchoLevel, 
                                              Vector& knot_vector_u_complete, Vector& knot_vector_v_complete, Vector& knot_vector_w_complete, double& knot_step_u, double& knot_step_v, double& knot_step_w,
                                              const Parameters refinements_parameters, const Parameters mParameters, ModelPart& surrogate_model_part_inner, 
                                              ModelPart& surrogate_model_part_outer);
                                              
        static void SnakeStep(ModelPart& skin_model_part, std::vector<std::vector<std::vector<int>>> &knot_spans_available, int idMatrixKnotSpansAvailable, 
                              int knot_span_u_1st_point, int knot_span_u_2nd_point, int knot_span_v_1st_point,int knot_span_v_2nd_point, 
                              double& x_true_boundary1, double& x_true_boundary2, double& y_true_boundary1, double& y_true_boundary2, 
                              double& knot_step_u, double& knot_step_v, const double starting_pos_u, const double starting_pos_v);
        
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

        static void MarkKnotSpansAvailable(std::vector<std::vector<std::vector<int>>>& knot_spans_available, int idInnerLoop, DynamicBins& testBins, 
                                           ModelPart& skin_model_part, double lambda, int insert_nb_per_span_u, int insert_nb_per_span_v, 
                                           double knot_step_u, double knot_step_v, const double starting_pos_u, const double starting_pos_v);

        static void CreateSurrogateBuondaryFromSnake_inner (std::vector<std::vector<std::vector<int>>> & knot_spans_available, int idInnerLoop, 
                            ModelPart& surrogate_model_part_inner, int insert_nb_per_span_u, int insert_nb_per_span_v, 
                            std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v,
                            const double starting_pos_u, const double starting_pos_v);

        static void CreateSurrogateBuondaryFromSnake_outer (DynamicBins& testBin_out, ModelPart& skin_model_part_out, 
                            std::vector<std::vector<std::vector<int>>>& knot_spans_available, int idMatrixKnotSpansAvailable , 
                            ModelPart& surrogate_model_part_outer, int insert_nb_per_span_u, int insert_nb_per_span_v, 
                            std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v,
                            const double starting_pos_u, const double starting_pos_v);


        static bool isPointInsideSkinBoundary3D(Point& pointToSearch, DynamicBins& testBins, ModelPart& skin_model_part_in);

        static void MarkKnotSpansAvailable3D( std::vector<std::vector<std::vector<std::vector<int>>>> & knot_spans_available, int idInnerLoop, DynamicBins& testBins, ModelPart& skin_model_part, double lambda, 
                                                int insert_nb_per_span_u, int insert_nb_per_span_v, int insert_nb_per_span_w, 
                                                double knot_step_u, double knot_step_v, double knot_step_w);

        static void CreateSurrogateBuondaryFromSnake_inner3D (std::vector<std::vector<std::vector<std::vector<int>>>> & knot_spans_available, int idInnerLoop, 
                                                                ModelPart& surrogate_model_part_inner, ModelPart& skin_model_part_in,  DynamicBins& testBins,
                                                                int insert_nb_per_span_u, int insert_nb_per_span_v, int insert_nb_per_span_w, 
                                                                std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v, std::vector<double>&  knot_vector_w,
                                                                double knot_step_u, double knot_step_v, double knot_step_w);

        static void CreateSurrogateBuondaryFromSnake_outer3D (DynamicBins& testBin_out, ModelPart& skin_model_part_out, std::vector<std::vector<std::vector<std::vector<int>>>>& knot_spans_available, 
                                                             int idMatrixKnotSpansAvailable , ModelPart& surrogate_model_part_outer,
                                                             int insert_nb_per_span_u, int insert_nb_per_span_v, int insert_nb_per_span_w, 
                                                             std::vector<double>& knot_vector_u, std::vector<double>&  knot_vector_v, std::vector<double>&  knot_vector_w);

    }; // Class DirectorUtilities

}  // namespace Kratos.

#endif // KRATOS_DIRECTOR_UTILITIES_H_INCLUDED  defined