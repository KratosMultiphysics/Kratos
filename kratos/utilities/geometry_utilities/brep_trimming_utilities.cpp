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
    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::CreateBrepSurfaceTrimmingIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops,
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        IntegrationInfo& rIntegrationInfo)
    {
        for (IndexType i_outer_loops = 0; i_outer_loops < rOuterLoops.size(); ++i_outer_loops) {

            Clipper2Lib::Paths64 all_loops(1 + rInnerLoops.size()), solution_outer, solution_inner;
            const double factor = 1e-10;

            Clipper2Lib::Point64 int_point;
            int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
            int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());
            for (IndexType j = 0; j < rOuterLoops[i_outer_loops].size(); ++j) {
                CurveTessellation<PointerVector<Node>> curve_tesselation;
                auto geometry_outer = *(rOuterLoops[i_outer_loops][j].get());
                curve_tesselation.Tessellate(
                    geometry_outer, 0.001, 1, true);
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

                int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
                int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());
                for (IndexType j = 0; j < rInnerLoops[i_inner_loops].size(); ++j) {
                    CurveTessellation<PointerVector<Node>> curve_tesselation;
                    auto geometry_inner = *(rInnerLoops[i_inner_loops][j].get());
                    curve_tesselation.Tessellate(
                        geometry_inner, 0.001, 1, true);
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
                    Clipper2Lib::Clipper64 clipper_operation_outer;
                    clipper_operation_outer.AddSubject(all_loops);

                    Clipper2Lib::Rect64 rectangle = Clipper2Lib::Rect64(
                        static_cast<cInt>(rSpansU[i] / factor), static_cast<cInt>(rSpansV[j] / factor),
                        static_cast<cInt>(rSpansU[i + 1] / factor), static_cast<cInt>(rSpansV[j + 1] / factor));

                    solution_outer = Clipper2Lib::RectClip(rectangle, all_loops);

                    const double span_area = std::abs(Clipper2Lib::Area(rectangle.AsPath()));
                    double clip_area = 0.0;
                    if (solution_outer.size() > 0)
                    {
                        clip_area = std::abs(Clipper2Lib::Area(solution_outer[0]));
                        for (IndexType k = 1; k < solution_outer.size(); ++k) {
                            clip_area -= std::abs(Clipper2Lib::Area(solution_outer[k]));
                        }
                    }

                    //operation for inner trimming
                    Clipper2Lib::Clipper64 clipper_operation_inner;
                    clipper_operation_inner.AddSubject(solution_outer);
                    clipper_operation_inner.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, solution_inner);

                    if (solution_outer.size() == 0 || clip_area/span_area < 1e-6) {
                        continue;
                    }
                    else if (std::abs(1-clip_area/span_area) < 1e-6) {
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
                    else {
                        std::vector<Matrix> triangles;
                        for(IndexType i = 0; i < solution_inner.size(); ++i)
                        {
                            BrepTrimmingUtilities::Triangulate_OPT(solution_inner[i], triangles, factor, span_area);
                        }

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
                    }
                    clipper_operation_outer.Clear();
                }
            }
        }
    };
    ///@} // Kratos Classes

    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::ComputeSpanTriangulation(
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops,
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops,
        const double u0,
        const double u1,
        const double v0,
        const double v1,
        bool& rIsTrimmed,
        std::vector<Matrix>& rTriangles)
    {
        using namespace Clipper2Lib;

        rTriangles.clear();
        rIsTrimmed = true; // assume trimmed unless proven otherwise

        const double factor = 1e-10;

        for (IndexType i_outer_loops = 0; i_outer_loops < rOuterLoops.size(); ++i_outer_loops)
        {
            Paths64 all_loops(1 + rInnerLoops.size());

            Point64 int_point;
            int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
            int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

            for (IndexType j = 0; j < rOuterLoops[i_outer_loops].size(); ++j)
            {
                CurveTessellation<PointerVector<Node>> curve_tesselation;
                auto geometry_outer = *(rOuterLoops[i_outer_loops][j].get());

                curve_tesselation.Tessellate(geometry_outer, 0.001, 1, true);
                auto tesselation = curve_tesselation.GetTessellation();

                for (IndexType u = 0; u < tesselation.size(); ++u)
                {
                    auto new_int_point = BrepTrimmingUtilities::ToIntPoint(
                        std::get<1>(tesselation[u])[0],
                        std::get<1>(tesselation[u])[1],
                        factor);

                    if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
                    {
                        all_loops[0].push_back(new_int_point); 
                        int_point = new_int_point;
                    }
                }
            }

            for (IndexType i_inner_loops = 0; i_inner_loops < rInnerLoops.size(); ++i_inner_loops)
            {
                int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
                int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

                for (IndexType j = 0; j < rInnerLoops[i_inner_loops].size(); ++j)
                {
                    CurveTessellation<PointerVector<Node>> curve_tesselation;
                    auto geometry_inner = *(rInnerLoops[i_inner_loops][j].get());

                    curve_tesselation.Tessellate(geometry_inner, 0.001, 1, true);
                    auto tesselation = curve_tesselation.GetTessellation();

                    for (IndexType u = 0; u < tesselation.size(); ++u)
                    {
                        auto new_int_point = BrepTrimmingUtilities::ToIntPoint(
                            std::get<1>(tesselation[u])[0],
                            std::get<1>(tesselation[u])[1],
                            factor);

                        if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
                        {
                            all_loops[i_inner_loops + 1].push_back(new_int_point);
                            int_point = new_int_point;
                        }
                    }
                }
            }

            // Clip with knot span (rectangle)
            Rect64 rectangle(
                static_cast<cInt>(u0 / factor),
                static_cast<cInt>(v0 / factor),
                static_cast<cInt>(u1 / factor),
                static_cast<cInt>(v1 / factor));

            Paths64 solution_outer = RectClip(rectangle, all_loops);

            const double span_area = std::abs(Area(rectangle.AsPath()));

            if (solution_outer.empty())
                continue;

            double clip_area = std::abs(Area(solution_outer[0]));
            for (IndexType k = 1; k < solution_outer.size(); ++k)
                clip_area -= std::abs(Area(solution_outer[k]));

            // Subtract inner loops
            Clipper64 clipper_operation_inner;
            clipper_operation_inner.AddSubject(solution_outer);

            Paths64 solution_inner;
            clipper_operation_inner.Execute(
                ClipType::Difference,
                FillRule::NonZero,
                solution_inner);

            // Classify the knot span
            if (clip_area / span_area < 1e-6)
            {
                continue; // empty for this outer loop
            }

            if (std::abs(1.0 - clip_area / span_area) < 1e-6)
            {
                // FULL → not trimmed
                rIsTrimmed = false;

                Matrix tri1(3,2), tri2(3,2);

                tri1(0,0)=u0; tri1(0,1)=v0;
                tri1(1,0)=u1; tri1(1,1)=v0;
                tri1(2,0)=u1; tri1(2,1)=v1;

                tri2(0,0)=u0; tri2(0,1)=v0;
                tri2(1,0)=u1; tri2(1,1)=v1;
                tri2(2,0)=u0; tri2(2,1)=v1;

                rTriangles.push_back(tri1);
                rTriangles.push_back(tri2);
            }
            else
            {
                rIsTrimmed = true;

                for (IndexType i = 0; i < solution_inner.size(); ++i)
                {
                    BrepTrimmingUtilities::Triangulate_OPT(
                        solution_inner[i],
                        rTriangles,
                        factor,
                        span_area);
                }
            }
        }
    }

    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::ComputeConstantTerms(
        Vector& rConstantTerms, IntegrationPointsArrayType& rElementIntegrationPoints,
        double U0, double U1, double V0, double V1,
        IntegrationInfo& rIntegrationInfo
        )
    {
        // Initialize const variables.
        IndexType ffactor = 1;
        const IndexType order_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1;
        const IndexType order_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1;

        // double check_number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor + 1);
        // if(rElementIntegrationPoints.size() < check_number_of_functions){
        //     ffactor = 1;
        // }

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
    template<bool TShiftedBoundary>
    double BrepTrimmingUtilities<TShiftedBoundary>::PointElimination(
        Vector& rConstantTerms, IntegrationPointsArrayType& rElementIntegrationPoints,
        IntegrationPointsArrayType& rElementNewIntegrationPoints,
        double U0, double U1, double V0, double V1,
        IntegrationInfo& rIntegrationInfo,
        const double clip_area
        )
    {
        /// Initialize variables.
        SizeType ffactor = 1;
        const SizeType order_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1;
        const SizeType order_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1;

        // double check_number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor + 1);
        // if(rElementIntegrationPoints.size() < check_number_of_functions){
        //     ffactor = 1;
        // }

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1);
        const IndexType min_number_of_points = order_u*order_v;

        const double targeted_residual = 1e-14; //1e-10, 1e-16 seg fault, to do
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
            global_residual = MomentFitting(rConstantTerms, rElementIntegrationPoints, U0, U1, V0, V1, rIntegrationInfo, clip_area);

            if( number_iterations == 0UL){
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

                // Remove points that are very small (< EPS1*max_value)
                // However, always keep #min_number_of_points.
                SizeType counter = 0;
                for(IndexType i = 0; i < rElementIntegrationPoints.size(); i++){
                    auto it = begin_it + i;
                    // TODO: Fix this > 2..4
                    if( it->Weight() < 1e-8*max_value && rElementIntegrationPoints.size() > min_number_of_points){
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

                // Leave loop in next iteration. Note if point_is_eliminated the moment fitting equation is solved again.
                if( rElementIntegrationPoints.size() <= min_number_of_points ){ //&& !point_is_eliminated ){
                    number_iterations = maximum_iteration + 10;
                }
            }
            number_iterations++;
        }
        if( (global_residual >= targeted_residual && prev_solution.size() > 0) || number_iterations > maximum_iteration ) {
            // Return previous solution.
            rElementNewIntegrationPoints.insert(rElementNewIntegrationPoints.begin(), prev_solution.begin(), prev_solution.end());
            rElementNewIntegrationPoints.erase(std::remove_if(rElementNewIntegrationPoints.begin(), rElementNewIntegrationPoints.end(), [](const IntegrationPointType& point) {
                return point.Weight() < 1e-14; }), rElementNewIntegrationPoints.end());

            return prev_residual;
        }
        else{
            // Return current solution.
            rElementNewIntegrationPoints.insert(rElementNewIntegrationPoints.begin(), rElementIntegrationPoints.begin(), rElementIntegrationPoints.end());
            rElementNewIntegrationPoints.erase(std::remove_if(rElementNewIntegrationPoints.begin(), rElementNewIntegrationPoints.end(), [](const IntegrationPointType& point) {
                return point.Weight() < 1e-14; }), rElementNewIntegrationPoints.end());

            return global_residual;
        }
    }
    template<bool TShiftedBoundary>
    double BrepTrimmingUtilities<TShiftedBoundary>::MomentFitting(
        Vector& rConstantTerms, IntegrationPointsArrayType& rElementIntegrationPoints,
        double U0, double U1, double V0, double V1,
        IntegrationInfo& rIntegrationInfo,
        const double clip_area
        )
    {
        // Initialize variables
        const double jacobian_x = 1.0;
        const double jacobian_y = 1.0;

        IndexType ffactor = 1;
        const IndexType order_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1;
        const IndexType order_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1;

        // double check_number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor + 1);
        // if(rElementIntegrationPoints.size() < check_number_of_functions){
        //     ffactor = 1;
        // }

        const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor + 1);
        const IndexType number_reduced_points = rElementIntegrationPoints.size();

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

        // Write computed weights onto integration points
        for( IndexType i = 0; i < number_reduced_points; ++i){
            // Divide by det_jacobian to account for the corresponding multiplication during the element integration within the used external solver.
            double new_weight = weights[i]/(jacobian_x*jacobian_y);
            rElementIntegrationPoints[i].SetWeight(new_weight);
        }

        return residual;    
    }

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::comp_heights_and_correction(std::list<c_vector<double,2> >& _polygon,
                                                 std::vector<c_vector<double,2> >& _eval_point_corr_range, std::vector<double>& _equal_height)
    {

        double Tol = 0.0001; //0.00001;

        std::list<c_vector<double,2> >::iterator start_iter = _polygon.begin();
        std::list<c_vector<double,2> >::iterator last_u_iter = _polygon.end();  //is one higher than the actual one
        last_u_iter--; //needs to be updated
        std::list<c_vector<double,2> >::iterator iter = _polygon.end();
        iter--;
        
        //container for saving polypoints for segment boundaries
        std::vector<std::list<c_vector<double,2> >::iterator> seg_bound_iterators;
        seg_bound_iterators.reserve(_equal_height.size()+1);
            
        //polygon should not be larger
        _eval_point_corr_range[0](1) = (*_polygon.begin())(1);
        seg_bound_iterators.push_back(_polygon.begin());

        bool is_bound_seg = false;
        bool last_found = false;
        double us = 0.0;

        //is coming from the left
        for(size_t i=1;i<_eval_point_corr_range.size()-1;i++) //loop over all evaluation points and segment points
        {
            if((i)%2==0)
            is_bound_seg=true;
            else
            is_bound_seg=false;

            double curren_eval_point =  _eval_point_corr_range[i](0);

            for (iter = start_iter; iter != last_u_iter; ++iter) 
            {
                iter++;
                if(iter == _polygon.end())
                {
                    return false;
                }
                    
                if((fabs((*iter)(0) +1.0)< Tol)&& !last_found)
                {
                    last_found = true;
                    iter++;
                    last_u_iter = iter;
                    iter--;
                }

                double next_poly_point = (*iter)(0);
                iter--;
                double current_poly_point = (*iter)(0);

                if((current_poly_point> curren_eval_point) && (curren_eval_point> next_poly_point))
                {
                    double current_h = (*iter)(1);
                    iter++;
                    double next_h = (*iter)(1);
                    //compute range in desired direction
                    double dv = next_h - current_h;
                    double du = next_poly_point - current_poly_point;
                    double u = curren_eval_point - current_poly_point;
                    c_vector<double,2> new_point;
                    new_point(0) = curren_eval_point;
                    new_point(1) = current_h + u/du * dv;
                    _eval_point_corr_range[i](1) = new_point(1);
                    _polygon.insert(iter,new_point);
                    iter--;  //to get to current
                    if(is_bound_seg)
                    seg_bound_iterators.push_back(iter);
                    start_iter = iter;
                    break;
                }
                else if ((next_poly_point == curren_eval_point))
                {
                    iter++;  //to next since tested against next
                    _eval_point_corr_range[i](1) = (*iter)(1);
                    if(is_bound_seg)
                    seg_bound_iterators.push_back(iter);
                    else
                    {
                    iter++;  
                    iter++;
                    last_u_iter = iter;    //test against next of the next run
                    iter--; //back to next
                    iter--;  //back to current
                    }
                    start_iter = iter;
                    break;
                }
            }
        }

        if(!last_found)
        {
            last_u_iter--;
            if(((*last_u_iter)(0) +1.0)< Tol)
            {
            last_u_iter++;
            }
            else
            return false;
        }

        last_u_iter--; //to get the correct one
        _eval_point_corr_range[_eval_point_corr_range.size()-1](1) = (*last_u_iter)(1);
        seg_bound_iterators.push_back(last_u_iter);
        
        if((_equal_height.size()+1) != seg_bound_iterators.size())
            return false;

        double ref_area = 0.0;
        for(size_t i=0;i<_equal_height.size();i++)
        {
            //compute delta s of segment
            double seg_u = fabs((*seg_bound_iterators[i])(0) - (*seg_bound_iterators[i+1])(0));
            //double incl = fabs((*seg_bound_iterators[i])(1) - (*seg_bound_iterators[i+1])(1))/(seg_u);
            double seg_area = 0.0;
            std::list<c_vector<double,2> >::iterator seg_end = seg_bound_iterators[i+1];
            std::list<c_vector<double,2> >::iterator seg_start = seg_bound_iterators[i];
            seg_end++;
            seg_start++;
            for(iter=seg_start;iter!=seg_end;iter++) //loop over all evaluation points and segment points
            {
            c_vector<double,2> current_poly = (*iter);
            iter--;
            c_vector<double,2> prev_poly = (*iter);
            iter++;
            double du = fabs(prev_poly(0) - current_poly(0));
            us += du;
            double dv = (prev_poly(1) + current_poly(1))/2.0 + 1.0;
            seg_area += du*dv;
            }
            ref_area += seg_area; 
        }

        double tol = 1e-12;
        double num_area=0.0;
        int counter = 0;
        std::vector<double> seg_us;
        std::vector<bool> buffer;
        buffer.resize(_equal_height.size()); //segments where error can be distributed
        seg_us.resize(_equal_height.size());
        for(size_t i=1;i<_eval_point_corr_range.size()-1;i=i+2)
        {
            buffer[counter] = true;
            //all evaluation points
            double segment_u = _eval_point_corr_range[i-1](0) - _eval_point_corr_range[i+1](0);
            double height = _eval_point_corr_range[i](1) + 1.0;
            if(height<=tol)
            {
            buffer[counter] = false;
            //height = 0.0;
            }
            else if (height>=2.0-tol)
            {
            buffer[counter] = false;
            //height = 2.0;
            }
            seg_us[counter] = segment_u;
            _equal_height[counter] = height;
            num_area+= height*segment_u;
            counter++;
        }

        double error = ref_area-num_area;

        //return true;
        int runs=0;
        for(runs=0;runs<4;runs++)
        {
            if(fabs(error)<tol)
                break;
            //bigger than zero and smaller than two
            counter=0;
            double length4error = 0.0;
            for(size_t i=0;i<_equal_height.size();i++)
            {
                //if(buffer[counter]==true)
                length4error += seg_us[counter];
                counter++;
            }
            counter=0;
            double dis_error =  error/length4error;
            //rInfo("mean_h %.10e",dis_error);

            //distribut error
            for(size_t i=0;i<_equal_height.size();i++)
            {
                //if(buffer[counter]==true)
                {
                _equal_height[counter] = _equal_height[counter] + dis_error;
                if(_equal_height[counter]<=tol)
                    {
                    //rInfo("hit_lower_boundary %.10e",_equal_height[counter]);
                    //buffer[counter] = false;
                    //_equal_height[counter] = 0.0;
                    }
                    else if (_equal_height[counter]>=2.0-tol)
                    {
                    //rInfo("hit_higer_boundary %.10e",_equal_height[counter]);
                    //buffer[counter] = false;
                    //_equal_height[counter] = 2.0;
                    }
                }
                counter++;
            }
            counter=0;
            //computer numerical error
            num_area=0.0;
            for(size_t i=0;i<_equal_height.size();i++)
            {
                num_area += _equal_height[counter]*seg_us[counter];
                counter++;
            }
            error = ref_area-num_area; //error must be distributed
        }
        return true;
    }


    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::check_polygon(std::list<c_vector<double,2> >& _polygon, c_vector<double,4> Borders, orientation& _replaced_border, 
                                                                std::list<c_vector<double,2> >::iterator& _iter)
    {   
        double Tol = 0.00001;

        _iter = _polygon.end();
        _iter--;
        c_vector<double,2> end_point = *_iter;
        c_vector<double,2> start_point = *_polygon.begin();

        ////Check whether the start and end point touch the borders, if not, this is not a valid curve
        orientation touching_start=od_none;
        orientation touching_end=od_none;

        double north = Borders(0);
        double east = Borders(1);
        double south = Borders(2);
        double west = Borders(3);
        double tolerance=0.00001;

        if(fabs(start_point(0)-west)<tolerance)
            touching_start=od_west;
        if(fabs(start_point(0)-east)<tolerance)
            touching_start=od_east;
        if(fabs(start_point(1)-north)<tolerance)
            touching_start=od_north;
        if(fabs(start_point(1)-south)<tolerance)
            touching_start=od_south;

        if(fabs(end_point(0)-west)<tolerance)
            touching_end=od_west;
        if(fabs(end_point(0)-east)<tolerance)
            touching_end=od_east;
        if(fabs(end_point(1)-north)<tolerance)
            touching_end=od_north;
        if(fabs(end_point(1)-south)<tolerance)
            touching_end=od_south;
        
        //trimming polygon is valid
        if(touching_start!=od_none && touching_end!=od_none) 
        {}
        else
            return false;

        //compute direction of trimming curve
        c_vector<double,2> g_trim_curve;
        g_trim_curve = end_point - start_point;

        //assign direction of curve to one of the boundaries such that inner part is on the left
        //enum orientation {od_none,od_north,od_east,od_south,od_west};
        _replaced_border = od_none; //edge which is replaced by polygon

        std::vector<c_vector<double,2> > edge_dir;
        edge_dir.resize(4);
        for(size_t i=0;i<edge_dir.size();i++)
            edge_dir[i].clear();

        //give directions
        edge_dir[0][0] = -1.0; //Borders[3] - Borders[1]; //north
        edge_dir[1][1] = 1.0; //Borders[0] - Borders[2];  //east
        edge_dir[2][0] = 1.0; //Borders[1] - Borders[3];  //south
        edge_dir[3][1] = -1.0; //Borders[2] - Borders[0]; //west
        

        double  best_value=0;
        std::vector<c_vector<double,2> > curves_dir; //main direction of curves
        std::vector<orientation> tmp_replace_edge;
        tmp_replace_edge.reserve(2);
        curves_dir.reserve(2);
        int best_index=-1;

        for(size_t i=0;i<4;i++) //determines the boundary which is replaced
        {
            double tmp_value = inner_prod(g_trim_curve,edge_dir[i]);
            if(tmp_value>0)
            {
            if(tmp_value>best_value)
            {
                best_index++;
                best_value = inner_prod(g_trim_curve,edge_dir[i])/norm_2(edge_dir[i]);
            }
            tmp_replace_edge.push_back(orientation((i+1)%5));
            curves_dir.push_back(edge_dir[i]);
            }
        }
        
        if(best_index==-1)
            return false;

        bool is_valid_curve = true;
        //check for all semgents of the trimmign curve
        for (_iter = _polygon.begin(); _iter != _polygon.end(); ++_iter) 
        {
            _iter++;
            if(_iter==_polygon.end())
            break;

            c_vector<double,2> seg_end = *_iter;
            _iter--;
            c_vector<double,2> seg_start = *_iter;

            c_vector<double,2> seg_tang = seg_end - seg_start;
            //computer tangent of current segment

            if(norm_2(seg_tang)<Tol)
            {
            _iter = _polygon.erase(_iter);
            continue;
            }

            double tmp = inner_prod(seg_tang,curves_dir[best_index]);
            if(tmp<-Tol)
            {
            is_valid_curve = false;
            break;
            }
        }

        //check other direction
        if(is_valid_curve)
            _replaced_border = tmp_replace_edge[best_index];
        else
        {
            is_valid_curve = true;
            best_index = (best_index+1)%2;
            for (_iter = _polygon.begin(); _iter != _polygon.end(); ++_iter) 
            {
                _iter++;
                if(_iter==_polygon.end())
                    break;

                c_vector<double,2> seg_end = *_iter;
                _iter--;
                c_vector<double,2> seg_start = *_iter;

                c_vector<double,2> seg_tang = seg_end - seg_start;
                //computer tangent of current segment

                if(norm_2(seg_tang)<Tol)
                {
                    _iter = _polygon.erase(_iter);
                    continue;
                }

                double tmp = inner_prod(seg_tang,curves_dir[best_index]);
                if(tmp<Tol)
                {
                    is_valid_curve = false;
                    break;
                }
            }
            if(is_valid_curve)
                _replaced_border = tmp_replace_edge[best_index];
        }
        return is_valid_curve;
    }

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::create_trimmed_domain(std::list<c_vector<double,2> >& _polygon, c_vector<double,4> Borders, orientation Replaced_border)
    {
        double tolerance=0.00001;
        //start point
        c_vector<double,2> poly_start = _polygon.front();
        c_vector<double,2> poly_end = _polygon.back();
        orientation touching_start=od_none;
        orientation touching_end=od_none;

        double north = Borders(0);
        double east = Borders(1);
        double south = Borders(2);
        double west = Borders(3);

            
        //  cornerpoints 
        //
        //            NORTH
        //         3 ------- 2
        //         |         |
        //   WEST  |         |  EAST
        //         |         |
        //         0 ------- 1
        //            SOUTH

        //check for start point touching
        if(fabs(poly_start(0)-west)<tolerance)
            touching_start=od_west;
        if(fabs(poly_start(0)-east)<tolerance)
            touching_start=od_east;
        if(fabs(poly_start(1)-north)<tolerance)
            touching_start=od_north;
        if(fabs(poly_start(1)-south)<tolerance)
            touching_start=od_south;

        if(fabs(poly_end(0)-west)<tolerance)
            touching_end=od_west;
        if(fabs(poly_end(0)-east)<tolerance)
            touching_end=od_east;
        if(fabs(poly_end(1)-north)<tolerance)
            touching_end=od_north;
        if(fabs(poly_end(1)-south)<tolerance)
            touching_end=od_south;

        if(touching_start!=od_none && touching_end!=od_none) //trimming polygon is valid
        {}
        else
            return false;

        //polygon needs to be extended
        std::vector<c_vector<double,2> > points;
        points.resize(4);
        //point 0
        points[0](0)=west;
        points[0](1)=south;
        //point 1
        points[1](0)=east;
        points[1](1)=south;
        //point 2
        points[2](0)=east;
        points[2](1)=north;
        //point 3
        points[3](0)=west;
        points[3](1)=north;

        
        //  cornerpoints 
        //
        //            NORTH
        //         3 ------- 2
        //         |         |
        //   WEST  |         |  EAST
        //         |         |
        //         0 ------- 1
        //            SOUTH

        if(Replaced_border==od_north) //
        { 
            if(touching_start==od_north)
            _polygon.push_front(points[2]);
            if(touching_end==od_north)
            _polygon.push_back(points[3]);
        }
        else if(Replaced_border==od_west)
        { 
            //push to upper or lower corner (left and right)
            if(touching_start==od_west)
            _polygon.push_front(points[3]);
            if(touching_end==od_west)
            _polygon.push_back(points[0]);
        }
        else if(Replaced_border==od_south)
        { 
            if(touching_start==od_south)
            _polygon.push_front(points[0]);
            if(touching_end==od_south)
            _polygon.push_back(points[1]);
        }
        else if(Replaced_border==od_east)
        { 
            if(touching_start==od_east)
            _polygon.push_front(points[1]);
            if(touching_end==od_east)
            _polygon.push_back(points[2]);
        }

        //start and end of polygon is on same side
        // if(touching_start==touching_end) //two cases
        // {
        //     return true;
        // }

        //search for possible next point at the end (in Gaussian space it starts always on the most left)
        int bound_vertex_index=0;
        if(touching_end==od_north)
            bound_vertex_index = 3;
        else if(touching_end==od_west)
            bound_vertex_index = 0;
        else if(touching_end==od_south)
            bound_vertex_index = 1;
        else if(touching_end==od_east)
            bound_vertex_index = 2;

        poly_start = _polygon.front();
        poly_end = _polygon.back();

        for(size_t i=0;i<4;i++)
        {
            int index = (i+bound_vertex_index)%4;

            if(is_on_Left(poly_start,poly_end,points[index]))
                _polygon.push_back(points[index]);
        }
        return true;
    }
    
    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::compute_bounding_box(std::list<c_vector<double,2> >& _polygon, c_vector<double,4>& _bounding_box)
    {
        std::list<c_vector<double,2> >::iterator iter = _polygon.begin();
        c_vector<double,2> tmp_point = *iter;

        _bounding_box[0] = tmp_point[1];
        _bounding_box[2] = tmp_point[1];
        _bounding_box[1] = tmp_point[0];
        _bounding_box[3] = tmp_point[0];

        //check for all semgents of the trimmign curve
        for (iter = _polygon.begin(); iter != _polygon.end(); ++iter) 
        {
            c_vector<double,2> tmp_point = *iter;

            if(tmp_point[1]>_bounding_box[0])//north
            _bounding_box[0] = tmp_point[1];

            if(tmp_point[0]>_bounding_box[1])//east
            _bounding_box[1] = tmp_point[0];

            if(tmp_point[1]<_bounding_box[2])//south
            _bounding_box[2] = tmp_point[1];
            
            if(tmp_point[0]<_bounding_box[3])//west
            _bounding_box[3] = tmp_point[0];
        }
    }

    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::map_polygon(std::list<c_vector<double,2> >& _polygon,
                                                 c_matrix<double,2,2> _rot, c_vector<double,2> _shifts, c_vector<double,2> _scales)
    {
        std::list<c_vector<double,2> >::iterator iter;
        for (iter = _polygon.begin(); iter != _polygon.end(); ++iter) 
        {
            c_vector<double,2> current_point = (*iter);
            current_point(0) = (current_point(0)-_shifts(0))*_scales(0);
            current_point(1) = (current_point(1)-_shifts(1))*_scales(1);
            axpy_prod(_rot, current_point, (*iter), true);
            //(*iter) = _rot * (*iter);
        }
    }

    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::map_quadrature_points(std::vector<std::vector<c_vector<double,2> > >& _quadpoints, 
                             c_matrix<double,2,2> _rot, c_vector<double,2> _shifts, c_vector<double,2> _scales)
    {
        std::list<c_vector<double,2> >::iterator iter;

        for(size_t i=0;i<_quadpoints.size();i++)
            for(size_t j=0;j<_quadpoints[i].size();j++)
            {
            c_vector<double,2> current_point;
            axpy_prod(_rot, _quadpoints[i][j],current_point, true);
            current_point(0) = (current_point(0)*_scales(0))+_shifts(0);
            current_point(1) = (current_point(1)*_scales(1))+_shifts(1);
            _quadpoints[i][j] = current_point;
            }
    }

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::is_on_Left(c_vector<double,2> a, c_vector<double,2> b, c_vector<double,2> c)
    {
        return ((b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0))) > 0;
    }

    template class KRATOS_API(KRATOS_CORE) BrepTrimmingUtilities<true>;
    template class KRATOS_API(KRATOS_CORE) BrepTrimmingUtilities<false>;
} // namespace Kratos.
