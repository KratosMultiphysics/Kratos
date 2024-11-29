//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// header includes
#include "brep_trimming_utilities.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    void BrepTrimmingUtilities::CreateBrepSurfaceTrimmingIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rOuterLoops,
        const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rInnerLoops,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        IntegrationInfo& rIntegrationInfo)
    {
        for (IndexType i_outer_loops = 0; i_outer_loops < rOuterLoops.size(); ++i_outer_loops) {

            Clipper2Lib::Paths64 all_loops(1 + rInnerLoops.size()), solution, solution_inner;
            const double factor = 1e-10;

            Clipper2Lib::Point64 int_point;
            int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
            int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());
            for (IndexType j = 0; j < rOuterLoops[i_outer_loops].size(); ++j) {
                CurveTessellation<PointerVector<Node>> curve_tesselation;
                auto geometry_outer = *(rOuterLoops[i_outer_loops][j].get());
                curve_tesselation.Tessellate(
                    geometry_outer, 1e-10, 1, true); ////// Tolerance of the polygon
                auto tesselation = curve_tesselation.GetTessellation();
                for (IndexType u = 0; u < tesselation.size(); ++u) {
                    auto new_int_point = BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y)) {
                        all_loops[i_outer_loops].push_back(new_int_point);
                        int_point.x = new_int_point.x;
                        int_point.y = new_int_point.y;
                    }
                }
            }

            for (IndexType i_inner_loops = 0; i_inner_loops < rInnerLoops.size(); ++i_inner_loops) {
                //ClipperLib::IntPoint int_point;
                int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
                int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());
                for (IndexType j = 0; j < rInnerLoops[i_inner_loops].size(); ++j) {
                    CurveTessellation<PointerVector<Node>> curve_tesselation;
                    auto geometry_inner = *(rInnerLoops[i_inner_loops][j].get());
                    curve_tesselation.Tessellate(
                        geometry_inner, 1e-10, 1, true);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (IndexType u = 0; u < tesselation.size(); ++u) {
                        auto new_int_point = BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                        if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y)) {
                            all_loops[i_inner_loops + 1].push_back(new_int_point);
                            int_point.x = new_int_point.x;
                            int_point.y = new_int_point.x;
                        }
                    }
                }
            }

            for (IndexType i = 0; i < rSpansU.size() - 1; ++i) {
                for (IndexType j = 0; j < rSpansV.size() - 1; ++j) {
                    Clipper2Lib::Clipper64 c;
                    c.AddSubject(all_loops);

                    Clipper2Lib::Rect64 rectangle = Clipper2Lib::Rect64(
                        static_cast<cInt>(rSpansU[i] / factor), static_cast<cInt>(rSpansV[j] / factor),
                        static_cast<cInt>(rSpansU[i + 1] / factor), static_cast<cInt>(rSpansV[j + 1] / factor));

                    solution = Clipper2Lib::RectClip(rectangle, all_loops);

                    const double span_area = std::abs(Clipper2Lib::Area(rectangle.AsPath()));
                    double clip_area = 0.0;  //if solution is empty, this crashes
                    if (solution.size() > 0)
                    {
                        clip_area = std::abs(Clipper2Lib::Area(solution[0]));
                        for (IndexType k = 1; k < solution.size(); ++k) {
                            clip_area -= std::abs(Clipper2Lib::Area(solution[k]));
                        }
                    }
                    // KRATOS_WATCH(clip_area)
                    
                    Clipper2Lib::Clipper64 d;
                    d.AddSubject(solution);// Add all input paths as subject polygons
                    // KRATOS_WATCH("I am outside clipper")
                    d.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, solution_inner);
                
                    if (solution.size() == 0 || clip_area == 0.0) {
                        continue;
                    }
                    else if (std::abs(clip_area- span_area) < 1000) {
                        const IndexType number_of_integration_points = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) * rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1);

                        IndexType initial_integration_size = rIntegrationPoints.size();

                        if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                            rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                        }

                        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        advance(integration_point_iterator, initial_integration_size);

                        IntegrationPointUtilities::IntegrationPoints2D(
                            integration_point_iterator,
                            rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1),
                            rSpansU[i], rSpansU[i + 1],
                            rSpansV[j], rSpansV[j + 1]);
                    }
                    else { // Here the actually cut the elements
                        std::vector<Matrix> triangles;
                        
                        // std::cout<<"i: "<<i<<", j: "<<j<<std::endl;
                        // KRATOS_WATCH(solution.size())
                        // for (IndexType i = 0; i < solution.size(); ++i)
                        // {   
                        //     for (size_t k = 0; k < solution[i].size(); ++k) {
                        //         std::cout << "Point solution " << i << "," << k << ": (" << BrepTrimmingUtilities::IntPointToDoublePoint(solution[i][k], factor) << ")" << std::endl;
                        //     }
                        // }
                        // reverse solution
                        // std::reverse(solution[1].begin(), solution[1].end());
                        // KRATOS_WATCH(std::abs(Clipper2Lib::Area(solution[1])))
                        // BrepTrimmingUtilities::Triangulate_OPT(solution[1], triangles, factor);
                   
                        for (IndexType i = 0; i < solution_inner.size(); ++i)
                        {   
                            KRATOS_WATCH(std::abs(Clipper2Lib::Area(solution_inner[i])))
                            BrepTrimmingUtilities::Triangulate_OPT(solution_inner[i], triangles, factor); 

                            
                                // // Print the points of the polygon
                                // for (size_t k = 0; k < solution_inner[i].size(); ++k) {
                                //     std::cout << "Point " << i << "," << k << ": (" << BrepTrimmingUtilities::IntPointToDoublePoint(solution_inner[i][k], factor) << ")" << std::endl;
                                // }
                        }
                        // KRATOS_WATCH(solution_inner.size())
                        // KRATOS_WATCH(triangles.size())

                        const SizeType number_of_points = std::max(rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1));

                        const IndexType number_of_integration_points = triangles.size() * IntegrationPointUtilities::s_gauss_triangle[number_of_points].size();

                        IndexType initial_integration_size = rIntegrationPoints.size();

                        if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                            rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                        }

                        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        advance(integration_point_iterator, initial_integration_size);

                        for (IndexType i = 0; i < triangles.size(); ++i)
                        {
                            IntegrationPointUtilities::IntegrationPointsTriangle2D(
                                integration_point_iterator,
                                number_of_points,
                                triangles[i](0, 0), triangles[i](1, 0), triangles[i](2, 0),
                                triangles[i](0, 1), triangles[i](1, 1), triangles[i](2, 1));
                        }
                        
                        // start moment fitting algorithm to reduce the number of integration points (modified from Meßmer et al, 2022)
                        // SizeType point_distribution_factor = 1.5; //gamma >= 2, TO DO
                        // const SizeType min_num_points = (rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0))*(rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1))*(point_distribution_factor);

                        // KRATOS_WATCH(rIntegrationPoints.size() - initial_integration_size)
                        // KRATOS_WATCH(min_num_points)
             
                        // if(((rIntegrationPoints.size() - initial_integration_size) >= min_num_points)) //do moment fitting (! to turn off)
                        // {
                        //     // std::cout<<"yuhu: "<<i<<", "<<j<<std::endl;

                        //     // Get integration points on the trimmed element.
                        //     IntegrationPointsArrayType element_integration_points;
                        //     IndexType element_integration_size = rIntegrationPoints.size() - initial_integration_size;
                        //     element_integration_points.resize(element_integration_size);

                        //     typename IntegrationPointsArrayType::iterator element_integration_point_iterator = element_integration_points.begin();
                        //     typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        //     advance(integration_point_iterator, initial_integration_size);

                        //     for (IndexType i = 0; i < element_integration_size; ++i)
                        //     {
                        //         (*element_integration_point_iterator) = (*integration_point_iterator);

                        //         element_integration_point_iterator++;
                        //         integration_point_iterator++;
                        //     }

                        //     // Get boundary integration points.
                        //     Vector constant_terms{};
                        //     ComputeConstantTerms(constant_terms, element_integration_points, rSpansU[i], rSpansU[i + 1], rSpansV[j], rSpansV[j + 1], rIntegrationInfo);

                        //     // KRATOS_WATCH(constant_terms)

                        //     // Start point elimination.
                        //     double residual = std::numeric_limits<double>::max();
                        //     IntegrationPointsArrayType element_new_integration_points;

                        //     // If residual can not be statisfied
                        //     if( residual > 1.0e-10){ 

                        //         // Run point elimination.
                        //         residual = PointElimination(constant_terms, element_integration_points, element_new_integration_points,
                        //                                     rSpansU[i], rSpansU[i + 1], rSpansV[j], rSpansV[j + 1], rIntegrationInfo);

                        //         // KRATOS_WATCH(residual)
                        //         // KRATOS_WATCH(element_new_integration_points.size())

                        //         // Erase the existing element integration points and replace it with the element new integration points
                        //         IndexType element_new_integration_size = element_new_integration_points.size();
                        //         // KRATOS_WATCH(rIntegrationPoints.size())
                        //         // KRATOS_WATCH(initial_integration_size)

                        //         rIntegrationPoints.erase(rIntegrationPoints.begin()+initial_integration_size, rIntegrationPoints.end());

                        //         // KRATOS_WATCH(rIntegrationPoints.size())
                        //         // KRATOS_WATCH(element_new_integration_size)

                        //         // rIntegrationPoints.resize(rIntegrationPoints.size() + element_new_integration_size);

                        //         rIntegrationPoints.insert(rIntegrationPoints.end(), element_new_integration_points.begin(), element_new_integration_points.end() );

                        //         // KRATOS_WATCH(rIntegrationPoints.size())
                        //     }
                        //     else
                        //     {
                        //         // TO DO
                        //     }
                        // }
                        // else
                        // {
                        //     continue;
                        // }
                    }
                    // c.Clear();
                }
            }
        }
    };

    void BrepTrimmingUtilities::ComputeConstantTerms(
        Vector& rConstantTerms, IntegrationPointsArrayType& rElementIntegrationPoints,
        double U0, double U1, double V0, double V1,
        IntegrationInfo& rIntegrationInfo
        )
    {
        // Initialize const variables.
        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1;
        const IndexType order_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1;

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor + 1);

        // Resize constant terms.
        rConstantTerms.resize(number_of_functions, false);
        std::fill( rConstantTerms.begin(),rConstantTerms.end(), 0.0);

        // Loop over all boundary integration points.
        IndexType row_index = 0UL;
        typename IntegrationPointsArrayType::iterator integration_point_iterator = rElementIntegrationPoints.begin();

        for (IndexType i = 0; i < rElementIntegrationPoints.size(); ++i)
        {
            // For all functions
            row_index = 0;
            const double weight = (*integration_point_iterator).Weight();
            for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
                for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                    // Assemble RHS
                    const double value = Polynomial::f_x((*integration_point_iterator)[0], i_x, U0, U1)
                                       * Polynomial::f_x((*integration_point_iterator)[1], i_y, V0, V1);
                    rConstantTerms[row_index] += value * weight;
                    row_index++;
                }
            }
            // Get next iterator
            integration_point_iterator++;
        }
    }

    double BrepTrimmingUtilities::PointElimination(
        Vector& rConstantTerms, IntegrationPointsArrayType& rElementIntegrationPoints,
        IntegrationPointsArrayType& rElementNewIntegrationPoints,
        double U0, double U1, double V0, double V1,
        IntegrationInfo& rIntegrationInfo
        )
    {
        /// Initialize variables.
        const SizeType ffactor = 1;
        const SizeType order_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1;
        const SizeType order_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1;
        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1);
        const IndexType min_number_of_points = order_u*order_v;

        const double targeted_residual = 1e-10;
        double global_residual = std::numeric_limits<double>::min();
        double prev_residual = 0.0;
        const SizeType maximum_iteration = 1000UL;
        SizeType number_iterations = 0UL;
        bool point_is_eliminated = false;
        IntegrationPointsArrayType prev_solution{};

        // If any point is eliminated, we must run another moment fitting loop, to guarantee that the weights are correct.
        // Also keep iterating, until targeted_residual is stepped over.
        while( point_is_eliminated || (global_residual < targeted_residual && number_iterations < maximum_iteration) ){
            point_is_eliminated = false;
            global_residual = MomentFitting(rConstantTerms, rElementIntegrationPoints, U0, U1, V0, V1, rIntegrationInfo);

            if( number_iterations == 0UL){
                // std::cout<<"first iteration"<<std::endl;
                /// In first iteration, revome all points but #number_of_functions
                // Sort integration points according to weight.
                std::sort(rElementIntegrationPoints.begin(), rElementIntegrationPoints.end(), [](const IntegrationPointType& point_a, const IntegrationPointType& point_b) -> bool {
                    return point_a.Weight() > point_b.Weight();
                });
                // Only keep #number_of_functions integration points.
                rElementIntegrationPoints.erase(rElementIntegrationPoints.begin()+number_of_functions, rElementIntegrationPoints.end());

                // Additionally remove all points that are zero.
                rElementIntegrationPoints.erase(std::remove_if(rElementIntegrationPoints.begin(), rElementIntegrationPoints.end(), [](const IntegrationPointType& point) {
                    return point.Weight() < 1e-14; }), rElementIntegrationPoints.end());

                // Stop if no points are left.
                if( rElementIntegrationPoints.size() == 0 )
                    break;

                point_is_eliminated = true;
            }
            else if( global_residual < targeted_residual && number_iterations < maximum_iteration){
                // std::cout<<"next iteration"<<std::endl;
                // Store solution, in previous solution
                prev_solution.clear();
                prev_solution.insert(prev_solution.begin(), rElementIntegrationPoints.begin(), rElementIntegrationPoints.end());
                prev_residual = global_residual;

                // Find min and max weights.
                auto min_value_it = rElementIntegrationPoints.begin();
                double min_value = std::numeric_limits<double>::max();
                double max_value = std::numeric_limits<double>::lowest();
                const auto begin_it = rElementIntegrationPoints.begin();
                for(IndexType i = 0; i < rElementIntegrationPoints.size(); i++){
                    auto it = begin_it + i;
                    if( it->Weight() < min_value ) {
                        min_value_it = it;
                        min_value = it->Weight();
                    }
                    if( it->Weight() > max_value ) {
                        max_value = it->Weight();
                    }
                }

                // KRATOS_WATCH(min_value)
                // KRATOS_WATCH(max_value)

                // Remove points that are very small (< EPS1*max_value)
                // However, always keep #min_number_of_points.
                SizeType counter = 0;
                for(IndexType i = 0; i < rElementIntegrationPoints.size(); i++){
                    auto it = begin_it + i;
                    // TODO: Fix this > 2..4
                    if( it->Weight() < 1e-8*max_value && rElementIntegrationPoints.size() > min_number_of_points){
                        // std::cout<<"remove i: "<<i<<std::endl;
                        rElementIntegrationPoints.erase(it);
                        point_is_eliminated = true;
                        counter++;
                    }
                }
                // If nothing was removed, remove at least one points.
                if( counter == 0 && rElementIntegrationPoints.size() > min_number_of_points){
                    rElementIntegrationPoints.erase(min_value_it);
                    point_is_eliminated = true;
                }

                // KRATOS_WATCH(rElementIntegrationPoints.size())

                // Leave loop in next iteration. Note if point_is_eliminated the moment fitting equation is solved again.
                if( rElementIntegrationPoints.size() <= min_number_of_points ){ //&& !point_is_eliminated ){
                    number_iterations = maximum_iteration + 10;
                }
            }
            number_iterations++;
        }
        if( (global_residual >= targeted_residual && prev_solution.size() > 0) || number_iterations > maximum_iteration ) {

            // std::cout << "prev_residual: " << prev_residual << std::endl;
            // Return previous solution.
            rElementNewIntegrationPoints.insert(rElementNewIntegrationPoints.begin(), prev_solution.begin(), prev_solution.end());
            rElementNewIntegrationPoints.erase(std::remove_if(rElementNewIntegrationPoints.begin(), rElementNewIntegrationPoints.end(), [](const IntegrationPointType& point) {
                return point.Weight() < 1e-14; }), rElementNewIntegrationPoints.end());

            // KRATOS_WATCH(rElementNewIntegrationPoints.size())

            return prev_residual;
        }
        else{
            // std::cout << "global_residual: " << global_residual << std::endl;
            // Return current solution.
            rElementNewIntegrationPoints.insert(rElementNewIntegrationPoints.begin(), rElementIntegrationPoints.begin(), rElementIntegrationPoints.end());
            rElementNewIntegrationPoints.erase(std::remove_if(rElementNewIntegrationPoints.begin(), rElementNewIntegrationPoints.end(), [](const IntegrationPointType& point) {
                return point.Weight() < 1e-14; }), rElementNewIntegrationPoints.end());

            // KRATOS_WATCH(rElementNewIntegrationPoints.size())

            return global_residual;
        }
    }

    double BrepTrimmingUtilities::MomentFitting(
        Vector& rConstantTerms, IntegrationPointsArrayType& rElementIntegrationPoints,
        double U0, double U1, double V0, double V1,
        IntegrationInfo& rIntegrationInfo
        )
    {
        // Initialize variables
        const double jacobian_x = 1.0;
        const double jacobian_y = 1.0;

        const IndexType ffactor = 1;
        const IndexType order_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1;
        const IndexType order_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1;

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor + 1);
        const IndexType number_reduced_points = rElementIntegrationPoints.size();

        // KRATOS_WATCH(number_reduced_points)

        /// Assemble moment fitting matrix.
        Matrix fitting_matrix(number_of_functions, number_reduced_points);
        IndexType row_index = 0;

        for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
            for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                // Loop over all points
                const auto points_it_begin = rElementIntegrationPoints.begin();
                for( IndexType column_index = 0; column_index < number_reduced_points; ++column_index ){
                    auto point_it = points_it_begin + column_index;

                    const double value = Polynomial::f_x((*point_it)[0], i_x, U0, U1)
                                       * Polynomial::f_x((*point_it)[1], i_y, V0, V1);

                    fitting_matrix(row_index,column_index) = value;
                }
                row_index++;
            }
        }

        // Solve non-negative Least-Square-Error problem.
        Vector weights(number_reduced_points);
        auto residual = nnls::nnls(fitting_matrix, rConstantTerms, weights)/number_of_functions;

        // KRATOS_WATCH(weights)

        // Write computed weights onto integration points
        for( IndexType i = 0; i < number_reduced_points; ++i){
            // Divide by det_jacobian to account for the corresponding multiplication during the element integration within the used external solver.
            double new_weight = weights[i]/(jacobian_x*jacobian_y);
            rElementIntegrationPoints[i].SetWeight(new_weight);
        }

        return residual;    
    }

    ///@} // Kratos Classes

    //template void BrepTrimmingUtilities::CreateBrepSurfaceTrimmingIntegrationPoints<
    //    DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>, Node>(
    //    BrepTrimmingUtilities::IntegrationPointsArrayType& rIntegrationPoints,
    //    const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rOuterLoops,
    //    const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rInnerLoops,
    //    const std::vector<double>& rSpansU,
    //    const std::vector<double>& rSpansV,
    //    IntegrationInfo& rIntegrationInfo);

} // namespace Kratos.