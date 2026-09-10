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
                        // Trimmed NURBS reparametrization 
                        std::vector<RuledSurfacePatch> patches;
                        const bool patch_ok = parametrize_local_trimmed_domain(
                            rOuterLoops, rInnerLoops, rSpansU[i], rSpansU[i + 1], rSpansV[j], rSpansV[j + 1], patches);

                        if (patch_ok && !patches.empty()) {
                            const SizeType npts_u = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0);
                            const SizeType npts_v = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1);

                            for (const auto& r_patch : patches) {
                                std::vector<double> patch_spans_u, patch_spans_v;
                                r_patch.Surface->SpansLocalSpace(patch_spans_u, 0);
                                r_patch.Surface->SpansLocalSpace(patch_spans_v, 1);

                                IntegrationPointsArrayType patch_integration_points(
                                    (patch_spans_u.size() - 1) * (patch_spans_v.size() - 1) * npts_u * npts_v);
                                auto patch_integration_point_iterator = patch_integration_points.begin();
                                for (IndexType su = 0; su + 1 < patch_spans_u.size(); ++su) {
                                    for (IndexType sv = 0; sv + 1 < patch_spans_v.size(); ++sv) {
                                        IntegrationPointUtilities::IntegrationPoints2D(
                                            patch_integration_point_iterator,
                                            npts_u, npts_v,
                                            patch_spans_u[su], patch_spans_u[su + 1],
                                            patch_spans_v[sv], patch_spans_v[sv + 1]);
                                    }
                                }

                                const IndexType initial_size = rIntegrationPoints.size();
                                rIntegrationPoints.resize(initial_size + patch_integration_points.size());

                                for (IndexType k = 0; k < patch_integration_points.size(); ++k) {
                                    const auto& r_local_point = patch_integration_points[k];

                                    array_1d<double, 2> position;
                                    double detJ;
                                    EvaluateRuledSurface(r_patch, r_local_point[0], r_local_point[1], position, detJ);

                                    rIntegrationPoints[initial_size + k][0] = position[0];
                                    rIntegrationPoints[initial_size + k][1] = position[1];
                                    rIntegrationPoints[initial_size + k].Weight() = r_local_point.Weight() * std::abs(detJ);
                                }
                            }
                        }

                        if (!patch_ok || patches.empty())
                        for(IndexType inner = 0; inner < solution_inner.size(); ++inner)
                        {
                            const auto& points_inner = solution_inner[inner];
                            auto on_border = [&rectangle](const Clipper2Lib::Point64& pt) {
                                return pt.x == rectangle.left || pt.x == rectangle.right ||
                                    pt.y == rectangle.top  || pt.y == rectangle.bottom;
                            };
                            IndexType excursions = 0;
                            bool prev_on_border = on_border(points_inner.back());  
                            for (const auto& pt : points_inner) {
                                bool cur_on_border = on_border(pt);
                                if (prev_on_border && !cur_on_border) ++excursions;  
                                prev_on_border = cur_on_border;
                            }
                            
                            bool is_agip = false;
                            orientation Replaced_border = od_none;
                            std::list<c_vector<double,2> > new_points;
                            c_vector<double,4> borders;
                            borders[0] = rectangle.bottom * factor;
                            borders[2] = rectangle.top * factor;
                            borders[1] = rectangle.right * factor;
                            borders[3] = rectangle.left * factor;

                            if(excursions == 1) //AGIP Candidates
                            {
                                //chech for validity of trimming curve if it can be integrated
                                std::list<c_vector<double,2> >::iterator iter;
                                std::vector< Clipper2Lib::Point64 > const& points = solution_inner[inner];
                                std::vector<Clipper2Lib::Point64> corners = {
                                    {rectangle.left, rectangle.top},{rectangle.left, rectangle.bottom},
                                    {rectangle.right, rectangle.top},{rectangle.right, rectangle.bottom}
                                };

                                std::vector<Clipper2Lib::Point64> filtered_points;
                                for (const auto& pt : points) {
                                    bool is_corner = false;
                                    for (const auto& corner : corners) {
                                        if (std::sqrt((pt.x - corner.x)*(pt.x - corner.x) + (pt.y - corner.y)*(pt.y - corner.y)) < 100) {
                                            is_corner = true;
                                            break;
                                        }
                                    }
                                    if (!is_corner) {
                                        filtered_points.push_back(pt);
                                    }
                                }

                                if(filtered_points.back() == points.back() && filtered_points.front() == points.front())
                                {
                                    Clipper2Lib::Point64 last = filtered_points.back();
                                    filtered_points.pop_back();
                                    filtered_points.insert(filtered_points.begin(), last);
                                }

                                for (const auto& pt : filtered_points) {
                                    c_vector<double, 2> vec;
                                    vec[0] = BrepTrimmingUtilities::IntPointToDoublePoint(pt, factor)[0];
                                    vec[1] = BrepTrimmingUtilities::IntPointToDoublePoint(pt, factor)[1];
                                    new_points.push_back(vec);
                                }

                                is_agip = check_polygon(new_points,borders,Replaced_border,iter);
                            }

                            if(is_agip == true)
                            {
                                create_trimmed_domain(new_points, borders, Replaced_border);
                            
                                std::vector<std::array<double, 2>> gauss_data_u = IntegrationPointUtilities::s_gauss_legendre[rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) - 1];
                                std::vector<std::array<double, 2>> gauss_data_v = IntegrationPointUtilities::s_gauss_legendre[rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) - 1];

                                std::vector<c_vector<double, 2>> gauss_u_w_weights;
                                gauss_u_w_weights.reserve(gauss_data_u.size());
                                for (const auto& point : gauss_data_u) {
                                    c_vector<double, 2> vec;
                                    vec[0] = point[0]*2-1;
                                    vec[1] = point[1]*2;
                                    gauss_u_w_weights.push_back(vec);
                                }
                                std::vector<c_vector<double, 2>> gauss_v_w_weights;
                                gauss_v_w_weights.reserve(gauss_data_v.size());
                                for (const auto& point : gauss_data_v) {
                                    c_vector<double, 2> vec;
                                    vec[0] = point[0]*2-1;
                                    vec[1] = point[1]*2;
                                    gauss_v_w_weights.push_back(vec);
                                }

                                if(Replaced_border==od_west || Replaced_border==od_east)
                                {
                                    std::vector<c_vector<double,2> > tmp = gauss_u_w_weights;
                                    gauss_u_w_weights = gauss_v_w_weights;
                                    gauss_v_w_weights = tmp;
                                }

                                c_vector<double,4> bounding_box;
                                compute_bounding_box(new_points,bounding_box);

                                double u1 = bounding_box(3); //west
                                double u2 = bounding_box(1); //east
                                double v1 = bounding_box(2); //north
                                double v2 = bounding_box(0); //south

                                //change borders
                                c_vector<double,4> Borders = bounding_box;

                                //mapping to Gauss quadrature space
                                c_vector<double,2> shifts;
                                shifts(0) = (u2 + u1)/2.0;
                                shifts(1) = (v2 + v1)/2.0;

                                c_vector<double,2> scals;
                                scals(0) = 2.0/(u2 - u1);
                                scals(1) = 2.0/(v2 - v1);

                                c_matrix<double,2,2> rot_matrix;
                                rot_matrix.clear();
                                if(Replaced_border==od_north)
                                {
                                    rot_matrix(0,0) = 1.0;
                                    rot_matrix(1,1) = 1.0;
                                }
                                else if(Replaced_border==od_east)
                                {
                                    rot_matrix(0,1) = -1.0;
                                    rot_matrix(1,0) = 1.0;
                                }
                                else if(Replaced_border==od_south)
                                {
                                    rot_matrix(0,0) = -1.0;
                                    rot_matrix(1,1) = -1.0;
                                }
                                else if(Replaced_border==od_west)
                                {
                                    rot_matrix(0,1) = 1.0;
                                    rot_matrix(1,0) = -1.0;
                                }
                                else
                                    break;
                                    // TODO
                                    // return quad_points;

                                //mapping in Gaussian space (shift, rotation, scaling)
                                map_polygon(new_points,rot_matrix,shifts,scals);

                                std::vector<double> eval_height(gauss_u_w_weights.size());
                                std::vector<c_vector<double,2> > u_eval_v_hight;
                                u_eval_v_hight.resize(2*gauss_u_w_weights.size()+1);

                                double tmp_start = 1.0;
                                u_eval_v_hight[0](0) = tmp_start;
                                int counter=gauss_u_w_weights.size()-1;
                                for(size_t i=1;i<u_eval_v_hight.size();i++)
                                {
                                    u_eval_v_hight[i](0) = gauss_u_w_weights[counter](0) + 1e-10; //Add 1e-10 to avoid numerical issues
                                    i++;
                                    u_eval_v_hight[i](0) = tmp_start - gauss_u_w_weights[counter](1); //Add 1e-10 to avoid numerical issues
                                    tmp_start = u_eval_v_hight[i](0);
                                    counter--;
                                }
                                u_eval_v_hight[u_eval_v_hight.size()-1](0) = -1.0;
                                
                                //Compute heights and correction
                                comp_heights_and_correction(new_points,u_eval_v_hight,eval_height);
                                // if(!comp_heights_and_correction(new_points,u_eval_v_hight,eval_height))
                                //     break;
                                //     // TODO :return quad_points;

                                std::vector<std::vector<c_vector<double,2> > > quadrature_points;
                                std::vector<std::vector<double> > int_weights;
                                quadrature_points.resize(gauss_u_w_weights.size());
                                int_weights.resize(gauss_u_w_weights.size());
                                for(size_t i=0;i<quadrature_points.size();i++)
                                {
                                    quadrature_points[i].resize(gauss_v_w_weights.size());
                                    int_weights[i].resize(gauss_v_w_weights.size());
                                }
                                //double sum=0.0;
                                counter = 0;
                                for(int i=quadrature_points.size()-1;i>=0;i--)
                                {
                                    //double sum_seg = 0.0;
                                    for(size_t j=0;j<quadrature_points[i].size();j++)
                                    {
                                    quadrature_points[i][j](0) = gauss_u_w_weights[i](0);
                                    double factor = eval_height[counter]*0.5; //length of g2 vector
                                    double tmp_height = (gauss_v_w_weights[j](0)+1.0)*eval_height[counter]*0.5;
                                    quadrature_points[i][j](1) = tmp_height-1.0;  //position is correct
                                    int_weights[i][j] = gauss_u_w_weights[i](1)*gauss_v_w_weights[j](1)*factor;

                                    }
                                    counter++;
                                }

                                //rotate quadrature points back
                                rot_matrix.clear();
                                if(Replaced_border==od_north)
                                {
                                    rot_matrix(0,0) = 1.0;
                                    rot_matrix(1,1) = 1.0;
                                }
                                else if(Replaced_border==od_east)
                                {
                                    rot_matrix(0,1) = 1.0;
                                    rot_matrix(1,0) = -1.0;
                                }
                                else if(Replaced_border==od_south)
                                {
                                    rot_matrix(0,0) = -1.0;
                                    rot_matrix(1,1) = -1.0;
                                }
                                else if(Replaced_border==od_west)
                                {
                                    rot_matrix(0,1) = -1.0;
                                    rot_matrix(1,0) = 1.0;
                                }
                                scals(0) = (u2 - u1)/2.0;
                                scals(1) = (v2 - v1)/2.0;

                                map_quadrature_points(quadrature_points,rot_matrix,shifts,scals);

                                double J3 = (u2-u1)*(v2-v1)*0.25;

                                // Assign the integration points
                                const IndexType number_of_integration_points = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) * rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1);

                                IndexType initial_integration_size = rIntegrationPoints.size();

                                if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                                    rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                                }

                                typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                                advance(integration_point_iterator, initial_integration_size);

                                for(size_t ku=0;ku<quadrature_points.size();ku++)
                                {
                                    for(size_t kv=0;kv<quadrature_points[ku].size();kv++)
                                    {
                                        (*integration_point_iterator)[0] = quadrature_points[ku][kv](0);
                                        (*integration_point_iterator)[1] = quadrature_points[ku][kv](1); 
                                        (*integration_point_iterator).Weight() = J3*int_weights[ku][kv];
                                        integration_point_iterator++;
                                    }
                                }
                            }
                            else //triangulation
                            {
                                std::vector<Matrix> triangles;
                                BrepTrimmingUtilities::Triangulate_OPT(solution_inner[inner], triangles, factor, span_area);

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
                                SizeType point_distribution_factor = 1.5; //1.5; //gamma >= 2, TO DO
                                const SizeType min_num_points = (rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0))*(rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1))*(point_distribution_factor);

                                bool do_moment_fitting = (rIntegrationPoints.size() - initial_integration_size) >= min_num_points;
                                if(do_moment_fitting)
                                {
                                    // Get integration points on the trimmed element.
                                    IntegrationPointsArrayType element_integration_points;
                                    IndexType element_integration_size = rIntegrationPoints.size() - initial_integration_size;
                                    element_integration_points.resize(element_integration_size);

                                    typename IntegrationPointsArrayType::iterator element_integration_point_iterator = element_integration_points.begin();
                                    typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                                    advance(integration_point_iterator, initial_integration_size);

                                    for (IndexType i = 0; i < element_integration_size; ++i)
                                    {
                                        (*element_integration_point_iterator) = (*integration_point_iterator);

                                        element_integration_point_iterator++;
                                        integration_point_iterator++;
                                    }

                                    // Get boundary integration points.
                                    Vector constant_terms{};
                                    ComputeConstantTerms(constant_terms, element_integration_points, rSpansU[i], rSpansU[i + 1], rSpansV[j], rSpansV[j + 1], rIntegrationInfo);

                                    // Start point elimination.
                                    double residual = std::numeric_limits<double>::max();
                                    IntegrationPointsArrayType element_new_integration_points;

                                    // If residual can not be statisfied
                                    if( residual > 1.0e-10){ 

                                        // Run point elimination.
                                        double clip_area_norm = clip_area * factor * factor;
                                        residual = PointElimination(constant_terms, element_integration_points, element_new_integration_points,
                                                                    rSpansU[i], rSpansU[i + 1], rSpansV[j], rSpansV[j + 1], rIntegrationInfo, clip_area_norm);

                                        // Erase the existing element integration points and replace it with the element new integration points
                                        IndexType element_new_integration_size = element_new_integration_points.size();

                                        rIntegrationPoints.erase(rIntegrationPoints.begin()+initial_integration_size, rIntegrationPoints.end());
                                        rIntegrationPoints.insert(rIntegrationPoints.end(), element_new_integration_points.begin(), element_new_integration_points.end());
                                    }
                                    else
                                    {
                                        KRATOS_WATCH("residual is small, no point elimination needed")
                                    }
                                }
                            }
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

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::CollectTrimCurveSegmentsForSpan(
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops,
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops,
        const double u0,
        const double u1,
        const double v0,
        const double v1,
        std::vector<TrimCurveSegment>& rTrimCurves)
    {
        rTrimCurves.clear();

        const double eps = 1e-9 * std::max(u1 - u0, v1 - v0);

        auto collect_from_loops = [&](const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rLoops)
        {
            for (IndexType i_loop = 0; i_loop < rLoops.size(); ++i_loop) {
                for (IndexType i_curve = 0; i_curve < rLoops[i_loop].size(); ++i_curve) {
                    const auto p_curve = rLoops[i_loop][i_curve];
                    const auto p_param_curve = p_curve->pGetCurveOnSurface()->pGetCurve();

                    // Breakpoints where this curve crosses the surface's knot-span grid.
                    std::vector<double> spans;
                    p_curve->SpansLocalSpace(spans);

                    for (IndexType k = 0; k + 1 < spans.size(); ++k) {
                        const double t_start = spans[k];
                        const double t_end = spans[k + 1];

                        CoordinatesArrayType local_mid = ZeroVector(3);
                        local_mid[0] = 0.5 * (t_start + t_end);
                        CoordinatesArrayType global_mid;
                        p_param_curve->GlobalCoordinates(global_mid, local_mid);

                        const bool inside_u = (global_mid[0] > u0 - eps) && (global_mid[0] < u1 + eps);
                        const bool inside_v = (global_mid[1] > v0 - eps) && (global_mid[1] < v1 + eps);

                        if (inside_u && inside_v) {
                            rTrimCurves.push_back({p_curve, t_start, t_end});
                        }
                    }
                }
            }
        };

        collect_from_loops(rOuterLoops);
        collect_from_loops(rInnerLoops);

        return true;
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::RawTrimCurve
    BrepTrimmingUtilities<TShiftedBoundary>::ExtractRawCurve(const TrimCurveSegment& rSegment)
    {
        const auto p_curve = rSegment.pCurve->pGetCurveOnSurface()->pGetCurve();

        RawTrimCurve raw;
        raw.Degree = p_curve->PolynomialDegree(0);
        raw.Knots = p_curve->Knots();
        raw.Weights = p_curve->Weights();
        raw.Points.resize(p_curve->size());
        for (IndexType i = 0; i < p_curve->size(); ++i) {
            const auto& point = (*p_curve)[i];
            raw.Points[i][0] = point[0];
            raw.Points[i][1] = point[1];
        }

        return RestrictCurveToRange(raw, rSegment.LocalParameterStart, rSegment.LocalParameterEnd);
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::RawTrimCurve
    BrepTrimmingUtilities<TShiftedBoundary>::RestrictCurveToRange(const RawTrimCurve& rCurve, const double TStart, const double TEnd)
    {
        const SizeType degree = rCurve.Degree;
        const double tol = 1e-9;

        auto count_multiplicity = [&](const Vector& rKnots, const double t) {
            SizeType mult = 0;
            for (IndexType i = 0; i < rKnots.size(); ++i) {
                if (std::abs(rKnots[i] - t) < tol) ++mult;
            }
            return mult;
        };

        const Vector full_knots = NurbsCurveRefinementUtilities::CreateFullKnotVector(rCurve.Knots, degree);

        std::vector<double> knots_to_insert;
        for (SizeType i = count_multiplicity(full_knots, TStart); i < degree + 1; ++i) knots_to_insert.push_back(TStart);
        for (SizeType i = count_multiplicity(full_knots, TEnd); i < degree + 1; ++i) knots_to_insert.push_back(TEnd);

        RawTrimCurve refined = rCurve;

        if (!knots_to_insert.empty()) {
            typename NurbsCurveGeometry<3, PointerVector<Node>>::PointsArrayType lifted_points;
            for (const auto& p : rCurve.Points) {
                lifted_points.push_back(Kratos::make_intrusive<Node>(0, p[0], p[1], 0.0));
            }

            const auto p_geometry_3d = (rCurve.Weights.size() > 0)
                ? Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(lifted_points, degree, rCurve.Knots, rCurve.Weights)
                : Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(lifted_points, degree, rCurve.Knots);

            PointerVector<Node> refined_points_3d;
            Vector refined_knots, refined_weights;
            NurbsCurveRefinementUtilities::KnotRefinement(*p_geometry_3d, knots_to_insert, refined_points_3d, refined_knots, refined_weights);

            refined.Knots = refined_knots;
            refined.Weights = (rCurve.Weights.size() > 0) ? refined_weights : Vector();
            refined.Points.resize(refined_points_3d.size());
            for (IndexType i = 0; i < refined_points_3d.size(); ++i) {
                refined.Points[i][0] = refined_points_3d(i)->X();
                refined.Points[i][1] = refined_points_3d(i)->Y();
            }
        }

        const Vector refined_full_knots = NurbsCurveRefinementUtilities::CreateFullKnotVector(refined.Knots, degree);

        auto first_index_of = [&](const double t) -> SizeType {
            for (IndexType i = 0; i < refined_full_knots.size(); ++i) {
                if (std::abs(refined_full_knots[i] - t) < tol) return i;
            }
            KRATOS_ERROR << "RestrictCurveToRange: parameter " << t << " not found in refined knot vector." << std::endl;
        };

        const SizeType cp_start = first_index_of(TStart);
        const SizeType k1 = first_index_of(TEnd);

        RawTrimCurve result;
        result.Degree = degree;
        result.Points.assign(refined.Points.begin() + cp_start, refined.Points.begin() + k1);

        if (refined.Weights.size() > 0) {
            result.Weights.resize(k1 - cp_start);
            for (IndexType i = 0; i < k1 - cp_start; ++i) result.Weights[i] = refined.Weights[cp_start + i];
        }

        Vector result_full_knots(k1 + degree + 1 - cp_start);
        for (IndexType i = 0; i < result_full_knots.size(); ++i) result_full_knots[i] = refined_full_knots[cp_start + i];
        result.Knots = NurbsCurveRefinementUtilities::CreateReducedKnotVector(result_full_knots, degree);

        return result;
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::RawTrimCurve
    BrepTrimmingUtilities<TShiftedBoundary>::DegreeElevateCurve(const RawTrimCurve& rCurve, const SizeType DegreeIncrease)
    {
        if (DegreeIncrease == 0) return rCurve;

        typename NurbsCurveGeometry<3, PointerVector<Node>>::PointsArrayType lifted_points;
        for (const auto& p : rCurve.Points) {
            lifted_points.push_back(Kratos::make_intrusive<Node>(0, p[0], p[1], 0.0));
        }

        const auto p_geometry_3d = (rCurve.Weights.size() > 0)
            ? Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(lifted_points, rCurve.Degree, rCurve.Knots, rCurve.Weights)
            : Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(lifted_points, rCurve.Degree, rCurve.Knots);

        PointerVector<Node> refined_points_3d;
        Vector refined_knots, refined_weights;
        NurbsCurveRefinementUtilities::DegreeElevation(*p_geometry_3d, DegreeIncrease, refined_points_3d, refined_knots, refined_weights);

        RawTrimCurve result;
        result.Degree = rCurve.Degree + DegreeIncrease;
        result.Knots = refined_knots;
        result.Weights = (rCurve.Weights.size() > 0) ? refined_weights : Vector();
        result.Points.resize(refined_points_3d.size());
        for (IndexType i = 0; i < refined_points_3d.size(); ++i) {
            result.Points[i][0] = refined_points_3d(i)->X();
            result.Points[i][1] = refined_points_3d(i)->Y();
        }
        return result;
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::RawTrimCurve
    BrepTrimmingUtilities<TShiftedBoundary>::ReverseCurve(const RawTrimCurve& rCurve)
    {
        RawTrimCurve reversed;
        reversed.Degree = rCurve.Degree;

        const SizeType n = rCurve.Points.size();
        reversed.Points.resize(n);
        for (IndexType i = 0; i < n; ++i) reversed.Points[i] = rCurve.Points[n - 1 - i];

        if (rCurve.Weights.size() > 0) {
            reversed.Weights.resize(n);
            for (IndexType i = 0; i < n; ++i) reversed.Weights[i] = rCurve.Weights[n - 1 - i];
        }

        const double t0 = rCurve.Knots[rCurve.Degree - 1];
        const double t1 = rCurve.Knots[rCurve.Knots.size() - rCurve.Degree];
        const SizeType m = rCurve.Knots.size();
        reversed.Knots.resize(m);
        for (IndexType i = 0; i < m; ++i) reversed.Knots[i] = t0 + t1 - rCurve.Knots[m - 1 - i];

        return reversed;
    }

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::MergeCurves(RawTrimCurve& rMasterCurve, RawTrimCurve CurveToAdd, const double Tolerance)
    {
        if (CurveToAdd.Degree > rMasterCurve.Degree) {
            rMasterCurve = DegreeElevateCurve(rMasterCurve, CurveToAdd.Degree - rMasterCurve.Degree);
        } else if (rMasterCurve.Degree > CurveToAdd.Degree) {
            CurveToAdd = DegreeElevateCurve(CurveToAdd, rMasterCurve.Degree - CurveToAdd.Degree);
        }
        const SizeType degree = rMasterCurve.Degree;

        auto close = [&](const array_1d<double, 2>& a, const array_1d<double, 2>& b) {
            return std::hypot(a[0] - b[0], a[1] - b[1]) < Tolerance;
        };

        if (close(rMasterCurve.Points.front(), CurveToAdd.Points.front()) ||
            close(rMasterCurve.Points.back(), CurveToAdd.Points.back())) {
            CurveToAdd = ReverseCurve(CurveToAdd);
        }

        const bool append_after = close(rMasterCurve.Points.back(), CurveToAdd.Points.front());
        const bool prepend_before = close(rMasterCurve.Points.front(), CurveToAdd.Points.back());

        if (!append_after && !prepend_before) return false;

        const Vector full_master = NurbsCurveRefinementUtilities::CreateFullKnotVector(rMasterCurve.Knots, degree);
        const Vector full_add = NurbsCurveRefinementUtilities::CreateFullKnotVector(CurveToAdd.Knots, degree);

        Vector weights_master;
        if (rMasterCurve.Weights.size() > 0) {
            weights_master = rMasterCurve.Weights;
        } else {
            weights_master.resize(rMasterCurve.Points.size());
            std::fill(weights_master.begin(), weights_master.end(), 1.0);
        }
        Vector weights_add;
        if (CurveToAdd.Weights.size() > 0) {
            weights_add = CurveToAdd.Weights;
        } else {
            weights_add.resize(CurveToAdd.Points.size());
            std::fill(weights_add.begin(), weights_add.end(), 1.0);
        }
        const bool is_rational = (rMasterCurve.Weights.size() > 0) || (CurveToAdd.Weights.size() > 0);

        RawTrimCurve merged;
        merged.Degree = degree;

        if (append_after) {
            const double last_knot = full_master[full_master.size() - 1];
            const double first_knot = full_add[0];

            merged.Points = rMasterCurve.Points;
            for (IndexType i = 1; i < CurveToAdd.Points.size(); ++i) merged.Points.push_back(CurveToAdd.Points[i]);

            std::vector<double> new_full;
            for (IndexType i = 0; i < full_master.size() - 1; ++i) new_full.push_back(full_master[i]);
            for (IndexType i = degree + 1; i < full_add.size(); ++i) new_full.push_back(full_add[i] - first_knot + last_knot);

            Vector new_full_v(new_full.size());
            for (IndexType i = 0; i < new_full.size(); ++i) new_full_v[i] = new_full[i];
            merged.Knots = NurbsCurveRefinementUtilities::CreateReducedKnotVector(new_full_v, degree);

            if (is_rational) {
                merged.Weights.resize(merged.Points.size());
                for (IndexType i = 0; i < weights_master.size(); ++i) merged.Weights[i] = weights_master[i];
                for (IndexType i = 1; i < weights_add.size(); ++i) merged.Weights[weights_master.size() - 1 + i] = weights_add[i];
            }
        } else {
            const double first_knot = full_master[0];
            const double last_knot = full_add[full_add.size() - 1];

            merged.Points = CurveToAdd.Points;
            for (IndexType i = 1; i < rMasterCurve.Points.size(); ++i) merged.Points.push_back(rMasterCurve.Points[i]);

            std::vector<double> new_full;
            for (IndexType i = 0; i < full_add.size() - 1; ++i) new_full.push_back(full_add[i]);
            for (IndexType i = degree + 1; i < full_master.size(); ++i) new_full.push_back(full_master[i] - first_knot + last_knot);

            Vector new_full_v(new_full.size());
            for (IndexType i = 0; i < new_full.size(); ++i) new_full_v[i] = new_full[i];
            merged.Knots = NurbsCurveRefinementUtilities::CreateReducedKnotVector(new_full_v, degree);

            if (is_rational) {
                merged.Weights.resize(merged.Points.size());
                for (IndexType i = 0; i < weights_add.size(); ++i) merged.Weights[i] = weights_add[i];
                for (IndexType i = 1; i < weights_master.size(); ++i) merged.Weights[weights_add.size() - 1 + i] = weights_master[i];
            }
        }

        rMasterCurve = merged;
        return true;
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::RawTrimCurve
    BrepTrimmingUtilities<TShiftedBoundary>::BuildStraightLineCurve(
        const array_1d<double, 2>& rPointA, const array_1d<double, 2>& rPointB,
        const double TStart, const double TEnd)
    {
        RawTrimCurve line;
        line.Degree = 1;
        line.Points = {rPointA, rPointB};
        line.Knots.resize(2);
        line.Knots[0] = TStart;
        line.Knots[1] = TEnd;
        return line;
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::RuledSurfacePatch
    BrepTrimmingUtilities<TShiftedBoundary>::BuildRuledSurfacePatch(
        const array_1d<double, 2>& rOppositeStart, const array_1d<double, 2>& rOppositeEnd,
        const RawTrimCurve& rTrimCurve)
    {
        const double t0 = rTrimCurve.Knots[rTrimCurve.Degree - 1];
        const double t1 = rTrimCurve.Knots[rTrimCurve.Knots.size() - rTrimCurve.Degree];

        RawTrimCurve opposite = BuildStraightLineCurve(rOppositeStart, rOppositeEnd, t0, t1);
        if (rTrimCurve.Degree > 1) {
            opposite = DegreeElevateCurve(opposite, rTrimCurve.Degree - 1);
        }

        const Vector full_trim = NurbsCurveRefinementUtilities::CreateFullKnotVector(rTrimCurve.Knots, rTrimCurve.Degree);
        std::vector<double> interior_knots;
        for (IndexType i = rTrimCurve.Degree + 1; i + rTrimCurve.Degree + 1 < full_trim.size(); ++i) {
            interior_knots.push_back(full_trim[i]);
        }

        if (!interior_knots.empty()) {
            typename NurbsCurveGeometry<3, PointerVector<Node>>::PointsArrayType lifted_points;
            for (const auto& p : opposite.Points) {
                lifted_points.push_back(Kratos::make_intrusive<Node>(0, p[0], p[1], 0.0));
            }
            const auto p_geometry_3d = Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(lifted_points, opposite.Degree, opposite.Knots);

            PointerVector<Node> refined_points_3d;
            Vector refined_knots, refined_weights;
            NurbsCurveRefinementUtilities::KnotRefinement(*p_geometry_3d, interior_knots, refined_points_3d, refined_knots, refined_weights);

            opposite.Knots = refined_knots;
            opposite.Points.resize(refined_points_3d.size());
            for (IndexType i = 0; i < refined_points_3d.size(); ++i) {
                opposite.Points[i][0] = refined_points_3d(i)->X();
                opposite.Points[i][1] = refined_points_3d(i)->Y();
            }
        }

        KRATOS_ERROR_IF(opposite.Points.size() != rTrimCurve.Points.size())
            << "BuildRuledSurfacePatch: opposite (" << opposite.Points.size()
            << ") and trim curve (" << rTrimCurve.Points.size() << ") control point counts do not match." << std::endl;

        const SizeType n = rTrimCurve.Points.size();
        typename ParametrizationPatchType::PointsArrayType surface_points;
        surface_points.resize(2 * n);
        for (IndexType v = 0; v < n; ++v) {
            surface_points(v * 2 + 0) = Kratos::make_shared<Point>(opposite.Points[v][0], opposite.Points[v][1]);
            surface_points(v * 2 + 1) = Kratos::make_shared<Point>(rTrimCurve.Points[v][0], rTrimCurve.Points[v][1]);
        }

        Vector knots_u(2);
        knots_u[0] = t0;
        knots_u[1] = t1;

        RuledSurfacePatch result;
        result.OppositeStart = rOppositeStart;
        result.OppositeEnd = rOppositeEnd;
        result.TrimCurve = rTrimCurve;

        if (rTrimCurve.Weights.size() > 0) {
            Vector weights(2 * n);
            for (IndexType v = 0; v < n; ++v) {
                weights[v * 2 + 0] = 1.0;
                weights[v * 2 + 1] = rTrimCurve.Weights[v];
            }
            result.Surface = Kratos::make_shared<ParametrizationPatchType>(surface_points, 1, rTrimCurve.Degree, knots_u, rTrimCurve.Knots, weights);
        } else {
            result.Surface = Kratos::make_shared<ParametrizationPatchType>(surface_points, 1, rTrimCurve.Degree, knots_u, rTrimCurve.Knots);
        }

        return result;
    }

    template<bool TShiftedBoundary>
    void BrepTrimmingUtilities<TShiftedBoundary>::EvaluateRuledSurface(
        const RuledSurfacePatch& rPatch,
        const double U, const double V,
        array_1d<double, 2>& rPosition,
        double& rDetJacobian)
    {
        const double t0 = rPatch.TrimCurve.Knots[rPatch.TrimCurve.Degree - 1];
        const double t1 = rPatch.TrimCurve.Knots[rPatch.TrimCurve.Knots.size() - rPatch.TrimCurve.Degree];
        const double length = t1 - t0;

        const double alpha = (V - t0) / length;
        array_1d<double, 2> opposite_pos;
        opposite_pos[0] = rPatch.OppositeStart[0] + alpha * (rPatch.OppositeEnd[0] - rPatch.OppositeStart[0]);
        opposite_pos[1] = rPatch.OppositeStart[1] + alpha * (rPatch.OppositeEnd[1] - rPatch.OppositeStart[1]);
        array_1d<double, 2> opposite_deriv;
        opposite_deriv[0] = (rPatch.OppositeEnd[0] - rPatch.OppositeStart[0]) / length;
        opposite_deriv[1] = (rPatch.OppositeEnd[1] - rPatch.OppositeStart[1]) / length;

        typename NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType trim_points;
        for (const auto& p : rPatch.TrimCurve.Points) {
            trim_points.push_back(Kratos::make_shared<Point>(p[0], p[1]));
        }
        NurbsCurveGeometry<2, PointerVector<Point>> trim_geometry = (rPatch.TrimCurve.Weights.size() > 0)
            ? NurbsCurveGeometry<2, PointerVector<Point>>(trim_points, rPatch.TrimCurve.Degree, rPatch.TrimCurve.Knots, rPatch.TrimCurve.Weights)
            : NurbsCurveGeometry<2, PointerVector<Point>>(trim_points, rPatch.TrimCurve.Degree, rPatch.TrimCurve.Knots);

        array_1d<double, 3> local_v; local_v[0] = V; local_v[1] = 0.0; local_v[2] = 0.0;
        std::vector<array_1d<double, 3>> trim_derivatives;
        trim_geometry.GlobalSpaceDerivatives(trim_derivatives, local_v, 1);

        const double N0 = (t1 - U) / length;
        const double N1 = (U - t0) / length;

        rPosition[0] = N0 * opposite_pos[0] + N1 * trim_derivatives[0][0];
        rPosition[1] = N0 * opposite_pos[1] + N1 * trim_derivatives[0][1];

        const double dXdU_x = (trim_derivatives[0][0] - opposite_pos[0]) / length;
        const double dXdU_y = (trim_derivatives[0][1] - opposite_pos[1]) / length;

        const double dXdV_x = N0 * opposite_deriv[0] + N1 * trim_derivatives[1][0];
        const double dXdV_y = N0 * opposite_deriv[1] + N1 * trim_derivatives[1][1];

        rDetJacobian = dXdU_x * dXdV_y - dXdU_y * dXdV_x;
    }

    template<bool TShiftedBoundary>
    typename BrepTrimmingUtilities<TShiftedBoundary>::Element_Face
    BrepTrimmingUtilities<TShiftedBoundary>::ClassifyFace(
        array_1d<double, 2>& rPoint,
        const double u0, const double u1, const double v0, const double v1,
        const double Tolerance)
    {
        const double diff_n = std::abs(rPoint[1] - v1);
        const double diff_e = std::abs(rPoint[0] - u1);
        const double diff_s = std::abs(rPoint[1] - v0);
        const double diff_w = std::abs(rPoint[0] - u0);

        Element_Face face = NONE;
        double best = 1e99;

        if (diff_n < Tolerance) { face = NORTH; rPoint[1] = v1; best = diff_n; }
        if (diff_e < Tolerance && diff_e < best) { face = EAST; rPoint[0] = u1; best = diff_e; }
        if (diff_s < Tolerance && diff_s < best) { face = SOUTH; rPoint[1] = v0; best = diff_s; }
        if (diff_w < Tolerance && diff_w < best) { face = WEST; rPoint[0] = u0; best = diff_w; }

        return face;
    }

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::calc_nurbs_patch(
        const std::vector<array_1d<double, 2>>& rCornerPoints,
        const double u0, const double u1, const double v0, const double v1,
        RawTrimCurve TrimCurve, Element_Face FaceStart, Element_Face FaceEnd,
        const double Tolerance,
        std::vector<RuledSurfacePatch>& rPatches)
    {
        const array_1d<double, 2>& SW = rCornerPoints[0];
        const array_1d<double, 2>& SE = rCornerPoints[1];
        const array_1d<double, 2>& NW = rCornerPoints[2];
        const array_1d<double, 2>& NE = rCornerPoints[3];

        const array_1d<double, 2> u_start_pt = TrimCurve.Points.front();
        const array_1d<double, 2> u_end_pt = TrimCurve.Points.back();
        const bool right_u_dir = (u_end_pt[0] - u_start_pt[0]) > 0;
        const bool right_v_dir = (u_end_pt[1] - u_start_pt[1]) > 0;

        array_1d<double, 2> opposite_a, opposite_b;

        // 2.1 curve intersects opposite faces
        if (FaceStart == SOUTH && FaceEnd == NORTH) { opposite_a = SW; opposite_b = NW; }
        else if (FaceStart == NORTH && FaceEnd == SOUTH) { opposite_a = NE; opposite_b = SE; }
        else if (FaceStart == WEST && FaceEnd == EAST) { opposite_a = NW; opposite_b = NE; }
        else if (FaceStart == EAST && FaceEnd == WEST) { opposite_a = SE; opposite_b = SW; }

        // 2.2 "triangular" patch -- clockwise adjacent
        else if (FaceStart == NORTH && FaceEnd == EAST) { opposite_a = NE; opposite_b = NE; }
        else if (FaceStart == EAST && FaceEnd == SOUTH) { opposite_a = SE; opposite_b = SE; }
        else if (FaceStart == SOUTH && FaceEnd == WEST) { opposite_a = SW; opposite_b = SW; }
        else if (FaceStart == WEST && FaceEnd == NORTH) { opposite_a = NW; opposite_b = NW; }

        // 2.3 "pentagon" patch -- counterclockwise adjacent, needs one straight extension
        else if (FaceStart == NORTH && FaceEnd == WEST) {
            RawTrimCurve extension = BuildStraightLineCurve(u_end_pt, SW, 0.0, 1.0);
            if (!MergeCurves(TrimCurve, extension, Tolerance)) return false;
            opposite_a = NE; opposite_b = SE;
        }
        else if (FaceStart == WEST && FaceEnd == SOUTH) {
            RawTrimCurve extension = BuildStraightLineCurve(u_end_pt, SE, 0.0, 1.0);
            if (!MergeCurves(TrimCurve, extension, Tolerance)) return false;
            opposite_a = NW; opposite_b = NE;
        }
        else if (FaceStart == SOUTH && FaceEnd == EAST) {
            RawTrimCurve extension = BuildStraightLineCurve(u_end_pt, NE, 0.0, 1.0);
            if (!MergeCurves(TrimCurve, extension, Tolerance)) return false;
            opposite_a = SW; opposite_b = NW;
        }
        else if (FaceStart == EAST && FaceEnd == NORTH) {
            RawTrimCurve extension = BuildStraightLineCurve(u_end_pt, NW, 0.0, 1.0);
            if (!MergeCurves(TrimCurve, extension, Tolerance)) return false;
            opposite_a = SE; opposite_b = SW;
        }

        // 2.4 curve intersects same face twice, needs two straight extensions
        else if (FaceStart == NORTH && FaceEnd == NORTH) {
            if (!right_u_dir) {
                RawTrimCurve front_ext = BuildStraightLineCurve(NE, u_start_pt, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, front_ext, Tolerance)) return false;
                RawTrimCurve back_ext = BuildStraightLineCurve(u_end_pt, NW, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, back_ext, Tolerance)) return false;
                opposite_a = SE; opposite_b = SW;
            } else {
                opposite_a = u_start_pt; opposite_b = u_end_pt;
            }
        }
        else if (FaceStart == EAST && FaceEnd == EAST) {
            if (right_v_dir) {
                RawTrimCurve front_ext = BuildStraightLineCurve(SE, u_start_pt, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, front_ext, Tolerance)) return false;
                RawTrimCurve back_ext = BuildStraightLineCurve(u_end_pt, NE, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, back_ext, Tolerance)) return false;
                opposite_a = SW; opposite_b = NW;
            } else {
                opposite_a = u_start_pt; opposite_b = u_end_pt;
            }
        }
        else if (FaceStart == SOUTH && FaceEnd == SOUTH) {
            if (right_u_dir) {
                RawTrimCurve front_ext = BuildStraightLineCurve(SW, u_start_pt, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, front_ext, Tolerance)) return false;
                RawTrimCurve back_ext = BuildStraightLineCurve(u_end_pt, SE, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, back_ext, Tolerance)) return false;
                opposite_a = NW; opposite_b = NE;
            } else {
                opposite_a = u_start_pt; opposite_b = u_end_pt;
            }
        }
        else if (FaceStart == WEST && FaceEnd == WEST) {
            if (!right_v_dir) {
                RawTrimCurve front_ext = BuildStraightLineCurve(NW, u_start_pt, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, front_ext, Tolerance)) return false;
                RawTrimCurve back_ext = BuildStraightLineCurve(u_end_pt, SW, 0.0, 1.0);
                if (!MergeCurves(TrimCurve, back_ext, Tolerance)) return false;
                opposite_a = NE; opposite_b = SE;
            } else {
                opposite_a = u_start_pt; opposite_b = u_end_pt;
            }
        }
        else {
            return false;
        }

        rPatches.push_back(BuildRuledSurfacePatch(opposite_a, opposite_b, TrimCurve));
        return true;
    }

    template<bool TShiftedBoundary>
    bool BrepTrimmingUtilities<TShiftedBoundary>::parametrize_local_trimmed_domain(
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops,
        const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops,
        const double u0,
        const double u1,
        const double v0,
        const double v1,
        std::vector<RuledSurfacePatch>& rPatches)
    {
        rPatches.clear();

        std::vector<TrimCurveSegment> trim_segments;
        if (!CollectTrimCurveSegmentsForSpan(rOuterLoops, rInnerLoops, u0, u1, v0, v1, trim_segments)) {
            return false;
        }

        if (trim_segments.empty()) {
            return true;
        }

        const double tol = 1e-9 * std::max(u1 - u0, v1 - v0);

        std::vector<RawTrimCurve> raw_curves;
        raw_curves.reserve(trim_segments.size());
        for (const auto& segment : trim_segments) {
            raw_curves.push_back(ExtractRawCurve(segment));
        }
        
        RawTrimCurve trim_curve = raw_curves[0];
        std::vector<bool> used(raw_curves.size(), false);
        used[0] = true;
        SizeType remaining = raw_curves.size() - 1;

        while (remaining > 0) {
            bool merged_any = false;
            for (IndexType k = 0; k < raw_curves.size(); ++k) {
                if (used[k]) continue;
                if (MergeCurves(trim_curve, raw_curves[k], tol)) {
                    used[k] = true;
                    --remaining;
                    merged_any = true;
                    break;
                }
            }
            if (!merged_any) {
                KRATOS_WATCH("Warning in parametrize_local_trimmed_domain: knot span is cut by more than one non-contiguous boundary arc. Cannot resolve a single trimming curve for this span.");
                return false;
            }
        }

        array_1d<double, 2> start_point = trim_curve.Points.front();
        array_1d<double, 2> end_point = trim_curve.Points.back();

        const Element_Face face_start = ClassifyFace(start_point, u0, u1, v0, v1, tol);
        const Element_Face face_end = ClassifyFace(end_point, u0, u1, v0, v1, tol);

        if (face_start == NONE || face_end == NONE) {
            return false;
        }

        trim_curve.Points.front() = start_point;
        trim_curve.Points.back() = end_point;

        std::vector<array_1d<double, 2>> corner_points(4);
        corner_points[0][0] = u0; corner_points[0][1] = v0; // SW
        corner_points[1][0] = u1; corner_points[1][1] = v0; // SE
        corner_points[2][0] = u0; corner_points[2][1] = v1; // NW
        corner_points[3][0] = u1; corner_points[3][1] = v1; // NE

        return calc_nurbs_patch(corner_points, u0, u1, v0, v1, trim_curve, face_start, face_end, tol, rPatches);
    }

    template class KRATOS_API(KRATOS_CORE) BrepTrimmingUtilities<true>;
    template class KRATOS_API(KRATOS_CORE) BrepTrimmingUtilities<false>;
} // namespace Kratos.
