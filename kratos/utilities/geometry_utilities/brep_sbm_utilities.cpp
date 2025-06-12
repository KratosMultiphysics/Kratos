//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// header includes
#include "brep_sbm_utilities.h"

namespace Kratos
{
    
///@name Kratos Classes
///@{
template<class TNodeType>
void BrepSbmUtilities<TNodeType>::CreateBrepSurfaceSbmIntegrationPoints(
    const std::vector<double>& rSpansU,
    const std::vector<double>& rSpansV,
    const GeometrySurrogateArrayType& rSurrogateOuterLoopGeometries,
    const GeometrySurrogateArrayType& rSurrogateInnerLoopGeometries,
    IntegrationPointsArrayType& rIntegrationPoints,
    IntegrationInfo& rIntegrationInfo)
{    
    // Loop over rSurrogateOuterLoopGeometries and collect the vertical conditions on each y-knot span
    std::vector<std::vector<double>> vertical_conditions_per_row(rSpansV.size()-1);
    const double tolerance = 1e-13;

    // Check if the outer loop is defined
    bool is_outer_loop_defined = false;
    if (rSurrogateOuterLoopGeometries.size() > 0) {is_outer_loop_defined = true;}
    
    // Loop over rSurrogateOuterLoopGeometries and collect the vertical conditions on each y-knot span
    for (auto i_cond : rSurrogateOuterLoopGeometries) {
        auto p_geometry = i_cond;
        if (std::abs((*p_geometry)[0].X()-(*p_geometry)[1].X()) < tolerance ) {
            // When local refining is performed, the surrogate steps might be larger that the refined knot spans
            const double direction = ((*p_geometry)[1].Y()-(*p_geometry)[0].Y())/std::abs(((*p_geometry)[1].Y()-(*p_geometry)[0].Y()));

            int i_row_1 = FindKnotSpans1D(rSpansV, (*p_geometry)[0].Y() + 1e-12*(direction));
            int i_row_2 = FindKnotSpans1D(rSpansV, (*p_geometry)[1].Y() - 1e-12*(direction));

            for (int i_row = std::min(i_row_1, i_row_2); i_row < std::max(i_row_1,i_row_2)+1; i_row++) {
                vertical_conditions_per_row[i_row].push_back((*p_geometry)[0].X());
            }
        }
    }
    
    // Loop over rSurrogateInnerLoopGeometries and collect the vertical conditions on each y-knot span
    for (auto i_cond : rSurrogateInnerLoopGeometries) {
        auto p_geometry = i_cond; 
        if (std::abs((*p_geometry)[0].X()-(*p_geometry)[1].X()) < tolerance ) {
            // When local refining is performed, the surrogate steps might be larger that the refined knot spans
            const double direction = ((*p_geometry)[1].Y()-(*p_geometry)[0].Y())/std::abs(((*p_geometry)[1].Y()-(*p_geometry)[0].Y()));
            int i_row_1 = FindKnotSpans1D(rSpansV, (*p_geometry)[0].Y() + 1e-12*(direction));
            int i_row_2 = FindKnotSpans1D(rSpansV, (*p_geometry)[1].Y() - 1e-12*(direction));

            for (int i_row = std::min(i_row_1, i_row_2); i_row < std::max(i_row_1,i_row_2)+1; i_row++) {
                vertical_conditions_per_row[i_row].push_back((*p_geometry)[0].X());
            }
        }
    }

    // Sort the vertical conditions by the x coordinate
    for (IndexType i_row = 0; i_row < rSpansV.size() - 1; ++i_row) {
        std::sort(vertical_conditions_per_row[i_row].begin(), vertical_conditions_per_row[i_row].end());  
    }   
    
    // Loop over each row and creation of the integration points
    for (IndexType j = 0; j < rSpansV.size() - 1; j++) {
        IndexType starting_column_index = 0;
        IndexType next_switch_knot_span;
        if (vertical_conditions_per_row[j].size() > 0) {
            next_switch_knot_span = FindKnotSpans1D(rSpansU, vertical_conditions_per_row[j][0]+tolerance);
            }
        else {
            next_switch_knot_span = rSpansU.size() - 1;
        }
        IndexType vertical_condition_index = 1;

        if (!is_outer_loop_defined) {
            vertical_condition_index = 1;
        } else {
            const auto& row = vertical_conditions_per_row[j];
            if (row.empty()) {
                continue;
            }
            starting_column_index = FindKnotSpans1D(rSpansU, row[0] + tolerance);
            next_switch_knot_span = FindKnotSpans1D(rSpansU, row[1] + tolerance);
            vertical_condition_index = 2;
        }
        
        for (IndexType i = starting_column_index; i < rSpansU.size() - 1; i++) {
            
            if (i == next_switch_knot_span) {
                if (vertical_condition_index < vertical_conditions_per_row[j].size()) {
                    i = FindKnotSpans1D(rSpansU, vertical_conditions_per_row[j][vertical_condition_index]+tolerance);
                } else {
                    break;
                }
                vertical_condition_index++;
                    if (vertical_condition_index < vertical_conditions_per_row[j].size()) {
                    next_switch_knot_span = FindKnotSpans1D(rSpansU, vertical_conditions_per_row[j][vertical_condition_index]+tolerance);
                }
                else {
                    next_switch_knot_span = rSpansU.size() - 1;
                }
                vertical_condition_index++;
            }            
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
        
    }
};

template<class TNodeType>
void BrepSbmUtilities<TNodeType>::CreateBrepVolumeSbmIntegrationPoints(
    IntegrationPointsArrayType& rIntegrationPoints,
    const std::vector<double>& rSpansU,
    const std::vector<double>& rSpansV,
    const std::vector<double>& rSpansW,
    GeometrySurrogateArrayType& rOuterLoops,
    GeometrySurrogateArrayType& rInnerLoops,
    IntegrationInfo& rIntegrationInfo)
{
    KRATOS_ERROR << "This function has not been implemented yet." << std::endl;
}

///@} // Kratos Classes
template<class TNodeType>
int BrepSbmUtilities<TNodeType>::FindKnotSpans1D(
    const std::vector<double>& rSpans, 
    const double coord) {
    
    for (IndexType i_row = 1; i_row < rSpans.size(); ++i_row) {
        if (rSpans[i_row] > coord) {
            return i_row-1;
        }
    }
    
    // If the coordinate lies at or beyond the last span boundary, return the last span index
    return rSpans.size()-1;
}

// TODO: 3D extension

template class KRATOS_API(KRATOS_CORE) BrepSbmUtilities<Node>;
} // namespace Kratos.
