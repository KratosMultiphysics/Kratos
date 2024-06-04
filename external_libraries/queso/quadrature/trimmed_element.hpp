//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef TRIMMED_ELEMENT_INCLUDE_HPP
#define TRIMMED_ELEMENT_INCLUDE_HPP

//// STL includes
#include <vector>
#include <array>
#include <variant>
#include <numeric>
//// Project includes
#include "queso/embedding/octree.h"
#include "queso/containers/element.hpp"
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/includes/parameters.h"
#include "queso/utilities/polynomial_utilities.h"
#include "queso/solvers/nnls.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  QuadratureTrimmedElement.
 * @author Manuel Messmer
 * @brief  Provides functions to create integration rules for trimmed elements.
 * @tparam TElementType
**/
template<typename TElementType>
class QuadratureTrimmedElement{
public:
    ///@name Type Definition
    ///@{
    typedef TElementType ElementType;
    typedef typename ElementType::IntegrationPointType IntegrationPointType;
    typedef typename ElementType::IntegrationPointVectorType IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef typename ElementType::BoundaryIntegrationPointType BoundaryIntegrationPointType;
    typedef std::vector<BoundaryIntegrationPointType> BoundaryIPsVectorType;
    typedef Unique<BoundaryIPsVectorType> BoundaryIPsVectorPtrType;
    typedef std::vector<double> VectorType;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Creates integration points for trimmed domain.
    ///@details 1. Distributes initial integration points uniformly in trimmed domain.
    ///         2. Computes constant terms of moment fitting equation.
    ///         3. Solves moment fitting equation in iterative point elimination algorithm.
    /// See: M. Meßmer et. al: Efficient CAD-integrated isogeometric analysis of trimmed solids,
    ///      Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
    ///@param rElement
    ///@param rIntegrationOrder
    ///@param Residual Targeted residual
    ///@param EchoLevel Default: 0
    static double AssembleIPs(ElementType& rElement, const Vector3i& rIntegrationOrder, double Residual, IndexType EchoLevel=0) {
        // Get boundary integration points.
        const auto p_trimmed_domain = rElement.pGetTrimmedDomain();
        const auto p_boundary_ips = p_trimmed_domain->template pGetBoundaryIps<typename TElementType::BoundaryIntegrationPointType>();

        // Get constant terms.
        VectorType constant_terms{};
        ComputeConstantTerms(constant_terms, p_boundary_ips, rElement, rIntegrationOrder);

        // Construct octree. Octree is used to distribute inital points within trimmed domain.
        const auto bounding_box = p_trimmed_domain->GetBoundingBoxOfTrimmedDomain();

        BoundingBoxType bounding_box_uvw = MakeBox( rElement.PointFromGlobalToParam(bounding_box.first),
                                                    rElement.PointFromGlobalToParam(bounding_box.second));

        Octree<TrimmedDomain> octree(p_trimmed_domain, bounding_box, bounding_box_uvw);

        // Start point elimination.
        double residual = MAXD;
        SizeType iteration = 0UL;
        SizeType point_distribution_factor = 1;
        IntegrationPointVectorType integration_points{};

        const IndexType max_iteration = (Math::Max(rIntegrationOrder) == 2) ? 4UL : 3UL;
        // If residual can not be statisfied, try with more points in initial set.
        while( residual > Residual && iteration < max_iteration){

            // Distribute intial points via an octree.
            const SizeType min_num_points = (rIntegrationOrder[0]+1)*(rIntegrationOrder[1]+1)*(rIntegrationOrder[2]+1)*(point_distribution_factor);
            DistributeIntegrationPoints(integration_points, octree, min_num_points, rIntegrationOrder);

            // If no point is contained in integration_points -> exit.
            if( integration_points.size() == 0 ){
                rElement.GetIntegrationPoints().clear();
                return 1;
            }

            // Also add old, moment fitted points to new set. 'old_integration_points' only contains points with weights > 0.0;
            auto& old_integration_points = rElement.GetIntegrationPoints();
            integration_points.insert(integration_points.end(), old_integration_points.begin(), old_integration_points.end() );
            old_integration_points.clear();

            // Run point elimination.
            residual = PointElimination(constant_terms, integration_points, rElement, rIntegrationOrder, Residual);

            // If residual is very high, remove all points. Note, elements without points will be neglected.
            if( residual > 1e-2 ) {
                auto& reduced_points = rElement.GetIntegrationPoints();
                reduced_points.clear();
            }

            // Update variables.
            point_distribution_factor *= 2;
            iteration++;
        }

        if( residual > Residual && EchoLevel > 2){
            QuESo_INFO << "Moment Fitting :: Targeted residual can not be achieved: " << residual << std::endl;
            // if( rParam.EchoLevel() > 3 ) {
            //     const std::string output_directory_name = rParam.Get<std::string>("output_directory_name");
            //     const std::string filename = output_directory_name + "/residual_not_achieved_id_" + std::to_string(rElement.GetId()) + ".stl";
            //     IO::WriteMeshToSTL(p_trimmed_domain->GetTriangleMesh(), filename.c_str(), true);
            // }
        }

        return residual;
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    /// @brief Distributes point within trimmed domain using an octree. In each leaf node, Gauss points according to rIntegrationOrder are generated.
    ///        Only points inside the trimmed domain are considered.
    ///        Every time this function is called the otree is refined and more points are distributed.
    /// @param[out] rIntegrationPoint
    /// @param rOctree
    /// @param MinNumPoints Minimum Number of Points
    /// @param rIntegrationOrder Order of Gauss quadrature.
    static void DistributeIntegrationPoints(IntegrationPointVectorType& rIntegrationPoint, Octree<TrimmedDomain>& rOctree, SizeType MinNumPoints, const Vector3i& rIntegrationOrder) {
        IndexType refinemen_level = rOctree.MaxRefinementLevel()+1;
        const IndexType max_iteration = 5UL;
        IndexType iteration = 0UL;
        while( rIntegrationPoint.size() < MinNumPoints && iteration < max_iteration){
            rOctree.Refine(std::min<IndexType>(refinemen_level, 4UL), refinemen_level);
            rIntegrationPoint.clear();
            rOctree.template AddIntegrationPoints<TElementType>(rIntegrationPoint, rIntegrationOrder);
            refinemen_level++;
            iteration++;
        }
    }

    /// @brief Computes constant terms of moment fitting equation via volume integration points.
    /// @param[out] rConstantTerms
    /// @param pIntegrationPoints (Unique<T>)
    /// @param rElement
    /// @param rIntegrationOrder
    static void ComputeConstantTerms(VectorType& rConstantTerms, const IntegrationPointVectorPtrType& pIntegrationPoints, const ElementType& rElement, const Vector3i& rIntegrationOrder) {
        // Initialize const variables.
        const PointType& a = rElement.GetBoundsUVW().first;
        const PointType& b = rElement.GetBoundsUVW().second;

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationOrder[0];
        const IndexType order_v = rIntegrationOrder[1];
        const IndexType order_w = rIntegrationOrder[2];

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);

        // Resize constant terms.
        rConstantTerms.resize(number_of_functions, false);
        std::fill( rConstantTerms.begin(),rConstantTerms.end(), 0.0);

        // Loop over all boundary integration points.
        IndexType row_index = 0UL;
        const auto begin_points_it = pIntegrationPoints->begin();
        for( IndexType i = 0; i < pIntegrationPoints->size(); ++i ){
            // Get iterator
            auto point_it = (begin_points_it + i);
            // For all functions
            row_index = 0;
            const double weight = point_it->Weight();
            for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
                for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                    for( IndexType i_z = 0; i_z <= order_w*ffactor; ++i_z){
                        // Assemble RHS
                        const double value = Polynomial::f_x(point_it->X(), i_x, a[0], b[0])
                            * Polynomial::f_x(point_it->Y(), i_y, a[1], b[1])
                            * Polynomial::f_x(point_it->Z(), i_z, a[2], b[2]);
                        rConstantTerms[row_index] += value * weight;
                        row_index++;
                    }
                }
            }
        }
    }

    /// @brief Computes constant terms of moment fitting equation via boundary integration points. This functions uses the divergence theorem
    //         to transform the respective volume integrals to countour/surface integrals.
    /// @param[out] rConstantTerms
    /// @param pBoundaryIPs (Unique<T>)
    /// @param rElement
    /// @param rIntegrationOrder
    static void ComputeConstantTerms(VectorType& rConstantTerms, const BoundaryIPsVectorPtrType& pBoundaryIPs, const ElementType& rElement, const Vector3i& rIntegrationOrder) {
        // Initialize const variables.
        const auto bounds_xyz = rElement.GetBoundsXYZ();

        // Constant terms / moments are evaluated in physical space.
        const PointType& a = bounds_xyz.first;
        const PointType& b = bounds_xyz.second;

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationOrder[0];
        const IndexType order_v = rIntegrationOrder[1];
        const IndexType order_w = rIntegrationOrder[2];

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);

        // Resize constant terms.
        rConstantTerms.resize(number_of_functions, false);
        std::fill( rConstantTerms.begin(),rConstantTerms.end(), 0.0);

        /// Initialize containers for f_x and f_x_int
        // X-direction
        std::vector<double> f_x_x(order_u*ffactor+1);
        std::vector<double> f_x_int_x(order_u*ffactor+1);
        // Y-direction
        std::vector<double> f_x_y(order_v*ffactor+1);
        std::vector<double> f_x_int_y(order_v*ffactor+1);
        // Z-direction
        std::vector<double> f_x_z(order_w*ffactor+1);
        std::vector<double> f_x_int_z(order_w*ffactor+1);

        // Loop over all boundary integration points.
        IndexType row_index = 0;
        const auto begin_points_it_ptr = pBoundaryIPs->begin();
        for( IndexType i = 0; i < pBoundaryIPs->size(); ++i ){
            // Note: The evaluation of polynomials is expensive. Therefore, precompute and store values
            // for f_x_x and f_x_int at each point.
            auto point_it = (begin_points_it_ptr + i);
            const auto& normal = point_it->Normal();
            PointType point = point_it->data();

            // X-Direction
            for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
                f_x_x[i_x] = Polynomial::f_x(point[0], i_x, a[0], b[0]);
                f_x_int_x[i_x] = Polynomial::f_x_int(point[0], i_x, a[0], b[0]);
            }
            // Y-Direction
            for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y){
                f_x_y[i_y] = Polynomial::f_x(point[1], i_y, a[1], b[1]);
                f_x_int_y[i_y] = Polynomial::f_x_int(point[1], i_y, a[1], b[1]);
            }
            // Z-Direction
            for( IndexType i_z = 0; i_z <= order_w*ffactor; ++i_z){
                f_x_z[i_z] = Polynomial::f_x(point[2], i_z, a[2], b[2]);
                f_x_int_z[i_z] = Polynomial::f_x_int(point[2], i_z, a[2], b[2]);
            }

            // Assembly RHS
            row_index = 0;
            const double weight = 1.0/3.0*point_it->Weight();
            for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
                for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                    for( IndexType i_z = 0; i_z <= order_w*ffactor; ++i_z){
                        // Compute normal for each face/triangle.
                        PointType value;
                        value[0] = f_x_int_x[i_x]*f_x_y[i_y]*f_x_z[i_z];
                        value[1] = f_x_x[i_x]*f_x_int_y[i_y]*f_x_z[i_z];
                        value[2] = f_x_x[i_x]*f_x_y[i_y]*f_x_int_z[i_z];

                        double integrand = normal[0]*value[0] + normal[1]*value[1] + normal[2]*value[2];
                        rConstantTerms[row_index] += integrand * weight;
                        row_index++;
                    }
                }
            }
        }
    }

    /// @brief Set-Up and solve moment fitting equation. Solve the moment fitting equation for the weights of the integration points.
    ///        Computed weights are directly assigned to rIntegrationPoint.
    /// @param rConstantTerms
    /// @param[out] rIntegrationPoint
    /// @param rElement
    /// @param rIntegrationOrder
    /// @return double Relative residual ||ax -b||_L2 / ||b||_L2
    static double MomentFitting(const VectorType& rConstantTerms, IntegrationPointVectorType& rIntegrationPoint, const ElementType& rElement, const Vector3i& rIntegrationOrder) {

        PointType a = rElement.GetBoundsUVW().first;
        PointType b = rElement.GetBoundsUVW().second;

        double jacobian = rElement.DetJ();

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationOrder[0];
        const IndexType order_v = rIntegrationOrder[1];
        const IndexType order_w = rIntegrationOrder[2];

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);
        const IndexType number_reduced_points = rIntegrationPoint.size();

        const double l2_norm_ct = std::sqrt( std::inner_product(rConstantTerms.begin(), rConstantTerms.end(), rConstantTerms.begin(), 0.0) );

        /// Assemble moment fitting matrix.
        NNLS::MatrixType fitting_matrix(number_of_functions * number_reduced_points, 0.0);
        const auto points_it_begin = rIntegrationPoint.begin();
        for( IndexType column_index = 0; column_index < number_reduced_points; ++column_index ){
            auto point_it = points_it_begin + column_index;
            IndexType row_index = 0;
            for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
                for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                    for( IndexType i_z = 0; i_z <= order_w*ffactor; ++i_z){
                        const double value = Polynomial::f_x(point_it->X(), i_x, a[0], b[0])
                                        * Polynomial::f_x(point_it->Y(), i_y, a[1], b[1])
                                        * Polynomial::f_x(point_it->Z(), i_z, a[2], b[2]);
                        // Matrix is serialized: Column first.
                        fitting_matrix[column_index*number_of_functions + row_index] = value;
                        row_index++;
                    }
                }
            }
        }

        // Solve non-negative Least-Square-Error problem.
        VectorType weights(number_reduced_points);
        VectorType tmp_constant_terms(rConstantTerms); // NNLS::solve does modify input. Therefore, copy is required.
        const double rel_residual = NNLS::solve(fitting_matrix, tmp_constant_terms, weights) / l2_norm_ct;

        // Write computed weights onto integration points
        for( IndexType i = 0; i < number_reduced_points; ++i){
            // Divide by det_jacobian to account for the corresponding multiplication during the element integration within the used external solver.
            double new_weight = weights[i]/(jacobian);
            rIntegrationPoint[i].SetWeight(new_weight);
        }

        return rel_residual;
    }

    /// @brief Start point elimination algorihtm. Final quadrature rule is stored in rElement.
    /// @param rConstantTerms
    /// @param rIntegrationPoint
    /// @param rElement
    /// @param rIntegrationOrder
    /// @param Residual targeted residual
    /// @return double achieved residual
    static double PointElimination(const VectorType& rConstantTerms, IntegrationPointVectorType& rIntegrationPoint, ElementType& rElement, const Vector3i& rIntegrationOrder, double Residual){
        /// Initialize variables.
        const SizeType ffactor = 1;
        const SizeType order_u = rIntegrationOrder[0];
        const SizeType order_v = rIntegrationOrder[1];
        const SizeType order_w = rIntegrationOrder[2];
        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);
        const IndexType min_number_of_points = order_u*order_v*order_w;

        const double targeted_residual = Residual;
        double global_residual = MIND;
        double prev_residual = 0.0;
        const SizeType maximum_iteration = 1000UL;
        SizeType number_iterations = 0UL;
        bool point_is_eliminated = false;
        IntegrationPointVectorType prev_solution{};

        // If any point is eliminated, we must run another moment fitting loop, to guarantee that the weights are correct.
        // Also keep iterating, until targeted_residual is stepped over.
        while( point_is_eliminated || (global_residual < targeted_residual && number_iterations < maximum_iteration) ){
            point_is_eliminated = false;
            global_residual = MomentFitting(rConstantTerms, rIntegrationPoint, rElement, rIntegrationOrder);
            if( number_iterations == 0UL){
                /// In first iteration, revome all points but #number_of_functions
                // Sort integration points according to weight.
                std::sort(rIntegrationPoint.begin(), rIntegrationPoint.end(), [](const IntegrationPointType& point_a, const IntegrationPointType& point_b) -> bool {
                        return point_a.Weight() > point_b.Weight();
                    });
                // Only keep #number_of_functions integration points.
                rIntegrationPoint.erase(rIntegrationPoint.begin()+number_of_functions, rIntegrationPoint.end());

                // Additionally remove all points that are zero.
                rIntegrationPoint.erase(std::remove_if(rIntegrationPoint.begin(), rIntegrationPoint.end(), [](const IntegrationPointType& point) {
                    return point.Weight() < ZEROTOL; }), rIntegrationPoint.end());

                // Stop if no points are left.
                if( rIntegrationPoint.size() == 0 )
                    break;

                point_is_eliminated = true;
            }
            else if( global_residual < targeted_residual && number_iterations < maximum_iteration){
                // Store solution, in previous solution
                prev_solution.clear();
                prev_solution.insert(prev_solution.begin(), rIntegrationPoint.begin(), rIntegrationPoint.end());
                prev_residual = global_residual;

                // Find min and max weights.
                auto min_value_it = rIntegrationPoint.begin();
                double min_value = MAXD;
                double max_value = LOWESTD;
                const auto begin_it = rIntegrationPoint.begin();
                for(IndexType i = 0; i < rIntegrationPoint.size(); i++){
                    auto it = begin_it + i;
                    if( it->Weight() < min_value ) {
                        min_value_it = it;
                        min_value = it->Weight();
                    }
                    if( it->Weight() > max_value ) {
                        max_value = it->Weight();
                    }
                }

                // Remove points that are very small (< EPS1*max_value)
                // However, always keep #min_number_of_points.
                SizeType counter = 0;
                for(IndexType i = 0; i < rIntegrationPoint.size(); i++){
                    auto it = begin_it + i;
                    // TODO: Fix this > 2..4
                    if( it->Weight() < 1e-8*max_value && rIntegrationPoint.size() > min_number_of_points){
                        rIntegrationPoint.erase(it);
                        point_is_eliminated = true;
                        counter++;
                    }
                }
                // If nothing was removed, remove at least one points.
                if( counter == 0 && rIntegrationPoint.size() > min_number_of_points){
                    rIntegrationPoint.erase(min_value_it);
                    point_is_eliminated = true;
                }
                // Leave loop in next iteration. Note if point_is_eliminated the moment fitting equation is solved again.
                if( rIntegrationPoint.size() <= min_number_of_points ){ //&& !point_is_eliminated ){
                    number_iterations = maximum_iteration + 10;
                }
            }
            number_iterations++;
        }

        auto& reduced_points = rElement.GetIntegrationPoints();
        if( (global_residual >= targeted_residual && prev_solution.size() > 0) || number_iterations > maximum_iteration ) {
            // Return previous solution.
            reduced_points.insert(reduced_points.begin(), prev_solution.begin(), prev_solution.end());
            reduced_points.erase(std::remove_if(reduced_points.begin(), reduced_points.end(), [](const IntegrationPointType& point) {
                return point.Weight() < ZEROTOL; }), reduced_points.end());

            return prev_residual;
        }
        else{
            // Return current solution.
            reduced_points.insert(reduced_points.begin(), rIntegrationPoint.begin(), rIntegrationPoint.end());
            reduced_points.erase(std::remove_if(reduced_points.begin(), reduced_points.end(), [](const IntegrationPointType& point) {
                return point.Weight() < ZEROTOL; }), reduced_points.end());

            return global_residual;
        }
    }
}; // End Class


} // End namespace queso

#endif // TRIMMED_ELEMENT_INCLUDE_HPP

// double QuadratureTrimmedElement::f_x(double x, int order, double a, double b){
//     double tmp_x = (2*x - a - b)/(b-a);
//     return p_n(tmp_x,order);
// }



// double QuadratureTrimmedElement::f_x(double x, int order){
//     if( order == 0){double
//     }
//     else {
//         return std::pow(x, order);
//     }
// }

// double QuadratureTrimmedElement::f_x_integral(double x, int order) {
//     if(  order == 0 ){
//         return x;
//     }
//     else {
//         return 1.0/ ( (double)order + 1.0) * std::pow(x, order+1);
//     }
// }

// double QuadratureTrimmedElement::f_x_integral(double x, int order, double a, double b){
//     switch(order)
//     {
//         case 0:
//             return x;
//         case 1:
//             return -std::pow((a + b - 2.0*x),2)/(4.0*(a - b));
//         case 2:
//             return - x/2.0 - std::pow((a + b - 2.0*x),3)/(4.0*std::pow((a - b),2));
//         case 3:
//             return (3.0*std::pow( (a + b - 2.0*x),2) )/(8.0*(a - b)) - (5*std::pow((a + b - 2.0*x),4))/(16*std::pow((a - b),3));
//         case 4:
//             return (3.0*x)/8.0 + (5.0*std::pow((a + b - 2.0*x),3))/(8*std::pow((a - b),2)) - (7.0*std::pow((a + b - 2*x),5))/(16*std::pow((a - b),4));
//         case 5:
//             return (35*std::pow((a + b - 2*x),4))/(32*std::pow((a - b),3)) - (15*std::pow((a + b - 2*x),2))/(32*(a - b)) - (21*std::pow((a + b - 2*x),6))/(32*std::pow((a - b),5));
//         case 6:
//             return (63*std::pow((a + b - 2*x),5))/(32*std::pow((a - b),4)) - (35*std::pow((a + b - 2*x),3))/(32*std::pow((a - b),2)) - (5*x)/16 - (33*std::pow((a + b - 2*x),7))/(32*std::pow((a - b),6));
//         case 7:
//             return (35.0*std::pow( (a + b - 2*x),2))/(64*(a - b)) - (315*std::pow((a + b - 2*x),4))/(128.0*std::pow((a - b),3)) + (231.0*std::pow((a + b - 2*x),6))/(64.0*std::pow((a - b),5)) - (429.0*std::pow((a + b - 2*x),8))/(256.0*std::pow((a - b),7));
//         case 8:
//             return (35.0*x)/128.0 + (105.0*std::pow( (a + b - 2*x),3) )/(64.0*std::pow( (a - b),2) ) - (693.0*std::pow( (a + b - 2*x),5) )/(128.0*std::pow((a - b),4) ) + (429.0*std::pow((a + b - 2*x),7))/(64.0*std::pow((a - b),6)) - (715.0*std::pow( (a + b - 2*x),9) )/(256.0*std::pow((a - b),8));
//     }

//     throw  std::invalid_argument("QuadratureTrimmedElement :: Order out of range!\n");
// }

// double QuadratureTrimmedElement::p_n(double x, int order) {
//     switch(order)
//     {
//         case 0:
//             return 1;
//         case 1:
//             return x;
//         case 2:
//             return 1.0/2.0*(3.0*std::pow(x,2)-1.0);
//         case 3:
//             return 1.0/2.0*(5.0*std::pow(x,3) - 3.0*x);
//         case 4:
//             return 1.0/8.0*(35.0*std::pow(x,4)-30.0*std::pow(x,2) +3.0);
//         case 5:
//             return 1.0/8.0*(63.0*std::pow(x,5)-70.0*std::pow(x,3)+15.0*x);
//         case 6:
//             return 1.0/16.0*(231.0*std::pow(x,6)-315.0*std::pow(x,4)+105.0*std::pow(x,2)-5.0);
//         case 7:
//             return 1.0/16.0*(429.0*std::pow(x,7)-693.0*std::pow(x,5)+315.0*std::pow(x,3)-35.0*x);
//         case 8:
//             return 1.0/128.0*(6435.0*std::pow(x,8) - 12012.0*std::pow(x,6)+6930.0*std::pow(x,4)-1260.0*std::pow(x,2)+35.0);
//     }

//     throw  std::invalid_argument("QuadratureTrimmedElement :: Order out of range!\n");
// }


