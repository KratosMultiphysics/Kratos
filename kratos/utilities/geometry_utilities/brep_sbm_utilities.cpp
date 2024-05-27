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

    void BrepSBMUtilities::CreateBrepSurfaceSBMIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        ModelPart& rSurrogateModelPart_inner, 
        ModelPart& rSurrogateModelPart_outer,
        IntegrationInfo& rIntegrationInfo)
    {
        // Loop over rSurrogateModelPart_outer and collect the vertical conditions on each y-knot span
        std::vector<std::vector<double>> verticalConditionsPerRow(rSpansV.size()-1);

        // Vector meshSizes_uv = rSurrogateModelPart_outer.GetProcessInfo().GetValue(MARKER_MESHES);
        bool isOuterLoopDefined = false;

        if (rSurrogateModelPart_outer.Conditions().size() > 0) {isOuterLoopDefined = true;}
        
        for (auto i_cond : rSurrogateModelPart_outer.Conditions()) {
            if (std::abs(i_cond.GetGeometry()[0].X()-i_cond.GetGeometry()[1].X()) < 1e-13 ) {
                //// Understand which knot span belongs to
                // const double mean_y_coord = (i_cond.GetGeometry()[0].Y() + i_cond.GetGeometry()[1].Y()) / 2;
                // int iRow = FindKnotSpans1D(rSpansV, mean_y_coord);
                //// Saving the x_coordinate of the current condition
                // verticalConditionsPerRow[iRow].push_back(i_cond.GetGeometry()[0].X());

                // When local refining is performed, the surrogate steps might be larger that the refined knot spans
                const double direction = (i_cond.GetGeometry()[1].Y()-i_cond.GetGeometry()[0].Y())/std::abs((i_cond.GetGeometry()[1].Y()-i_cond.GetGeometry()[0].Y()));

                int iRow1 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[0].Y() + 1e-12*(direction));
                int iRow2 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[1].Y() - 1e-12*(direction));

                for (int i_row = std::min(iRow1, iRow2); i_row < std::max(iRow1,iRow2)+1; i_row++) {
                    verticalConditionsPerRow[i_row].push_back(i_cond.GetGeometry()[0].X());
                }
            }
        }
        
        for (auto i_cond : rSurrogateModelPart_inner.Conditions()) {
            if (std::abs(i_cond.GetGeometry()[0].X()-i_cond.GetGeometry()[1].X()) < 1e-13 ) {
                // // Understand which knot span belongs to
                // const double mean_y_coord = (i_cond.GetGeometry()[0].Y() + i_cond.GetGeometry()[1].Y()) / 2;
                // int iRow = FindKnotSpans1D(rSpansV, mean_y_coord);
                // // Saving the x_coordinate of the current condition
                // verticalConditionsPerRow[iRow].push_back(i_cond.GetGeometry()[0].X());

                // When local refining is performed, the surrogate steps might be larger that the refined knot spans
                const double direction = (i_cond.GetGeometry()[1].Y()-i_cond.GetGeometry()[0].Y())/std::abs((i_cond.GetGeometry()[1].Y()-i_cond.GetGeometry()[0].Y()));
                int iRow1 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[0].Y() + 1e-12*(direction));
                int iRow2 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[1].Y() - 1e-12*(direction));

                for (int i_row = std::min(iRow1, iRow2); i_row < std::max(iRow1,iRow2)+1; i_row++) {
                    verticalConditionsPerRow[i_row].push_back(i_cond.GetGeometry()[0].X());
                }
            }
        }

        // Sort the vertical conditions by the x coordinate
        for (IndexType iRow = 0; iRow < rSpansV.size() - 1; ++iRow) {
            // verticalConditionsPerRow[iRow].push_back(0.0);
            // verticalConditionsPerRow[iRow].push_back(2.0);
            std::sort(verticalConditionsPerRow[iRow].begin(), verticalConditionsPerRow[iRow].end());  
            // Remove duplicates
            // auto last = std::unique(verticalConditionsPerRow[iRow].begin(), verticalConditionsPerRow[iRow].end());
            // verticalConditionsPerRow[iRow].erase(last, verticalConditionsPerRow[iRow].end());
        }   
        // KRATOS_WATCH(verticalConditionsPerRow)
        
        // Loop over each row and creation of the integration points
        for (IndexType j = 0; j < rSpansV.size() - 1; j++) {
            int startingColumnIndex = 0;
            int nextSwitchKnotSpan;
            if (verticalConditionsPerRow[j].size() > 0) {
                // nextSwitchKnotSpan = (verticalConditionsPerRow[j][0]+1e-13)/meshSizes_uv[0];
                nextSwitchKnotSpan = FindKnotSpans1D(rSpansU, verticalConditionsPerRow[j][0]+1e-13);
                }
            else {
                nextSwitchKnotSpan = rSpansU.size() - 1;
            }
            int verticalConditionIndex = 1;

            if (isOuterLoopDefined){
                if (verticalConditionsPerRow[j].empty()) {continue;}
                // startingColumnIndex = (verticalConditionsPerRow[j][0]+1e-13)/meshSizes_uv[0];
                startingColumnIndex = FindKnotSpans1D(rSpansU, verticalConditionsPerRow[j][0]+1e-13);
                // nextSwitchKnotSpan = (verticalConditionsPerRow[j][1]+1e-13)/meshSizes_uv[0];
                nextSwitchKnotSpan = FindKnotSpans1D(rSpansU, verticalConditionsPerRow[j][1]+1e-13);
                verticalConditionIndex = 2;
                }
            for (IndexType i = startingColumnIndex; i < rSpansU.size() - 1; i++) {
                
                if (i == nextSwitchKnotSpan) {
                    if (verticalConditionIndex < verticalConditionsPerRow[j].size()) {
                        i = FindKnotSpans1D(rSpansU, verticalConditionsPerRow[j][verticalConditionIndex]+1e-13);
                    } else {
                        break;
                    }
                    verticalConditionIndex++;
                     if (verticalConditionIndex < verticalConditionsPerRow[j].size()) {
                        nextSwitchKnotSpan = FindKnotSpans1D(rSpansU, verticalConditionsPerRow[j][verticalConditionIndex]+1e-13);
                    }
                    else {
                        nextSwitchKnotSpan = rSpansU.size() - 1;
                    }
                    verticalConditionIndex++;
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
    ///@} // Kratos Classes

    int BrepSBMUtilities::FindKnotSpans1D(const std::vector<double>& rSpans, const double coord) {
        
        for (IndexType i_row = 1; i_row < rSpans.size(); ++i_row) {
            if (rSpans[i_row] > coord) {
                return i_row-1;
            }
        }
    }

    // 3D
    void BrepSBMUtilities::CreateBrepVolumeSBMIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        const std::vector<double>& rSpansW,
        ModelPart& rSurrogateModelPart_inner, 
        ModelPart& rSurrogateModelPart_outer,
        IntegrationInfo& rIntegrationInfo)
    {
        // Loop over rSurrogateModelPart_inner and outer and save the perpendicular faces with respect to x-direction
        std::vector<std::vector<std::vector<double>>> perpendicular_conditions_per_u_direction(rSpansV.size()-1);

        for (int v = 0; v < rSpansV.size()-1; v++) 
        {
            perpendicular_conditions_per_u_direction[v].resize(rSpansV.size()-1);
        }
        bool isOuterLoopDefined = false;
        if (rSurrogateModelPart_outer.Conditions().size() > 0) {isOuterLoopDefined = true;}
        
        for (auto i_cond : rSurrogateModelPart_outer.Conditions()) {
            // Check if the condition is perpendicular to x-direction
            if (std::abs(i_cond.GetGeometry()[0].X()-i_cond.GetGeometry()[2].X()) < 1e-13 ) {

                // When local refining is performed, the surrogate steps might be larger that the refined knot spans
                const double direction_v = (i_cond.GetGeometry()[2].Y()-i_cond.GetGeometry()[0].Y())/std::abs((i_cond.GetGeometry()[2].Y()-i_cond.GetGeometry()[0].Y()));
                int vRow1 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[0].Y() + 1e-12*(direction_v));
                int vRow2 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[2].Y() - 1e-12*(direction_v));

                const double direction_w = (i_cond.GetGeometry()[2].Z()-i_cond.GetGeometry()[0].Z())/std::abs((i_cond.GetGeometry()[2].Z()-i_cond.GetGeometry()[0].Z()));
                int wRow1 = FindKnotSpans1D(rSpansW, i_cond.GetGeometry()[0].Z() + 1e-12*(direction_w));
                int wRow2 = FindKnotSpans1D(rSpansW, i_cond.GetGeometry()[2].Z() - 1e-12*(direction_w));

                for (int v_row = std::min(vRow1, vRow2); v_row < std::max(vRow1,vRow2)+1; v_row++) {
                    for (int w_row = std::min(wRow1, wRow2); w_row < std::max(wRow1,wRow2)+1; w_row++) {
                        perpendicular_conditions_per_u_direction[v_row][w_row].push_back(i_cond.GetGeometry()[0].X());
                    }
                }
            }
        }
        
        for (auto i_cond : rSurrogateModelPart_inner.Conditions()) {
            // Check if the condition is perpendicular to x-direction
            if (std::abs(i_cond.GetGeometry()[0].X()-i_cond.GetGeometry()[2].X()) < 1e-13 ) {

                // When local refining is performed, the surrogate steps might be larger that the refined knot spans
                const double direction_v = (i_cond.GetGeometry()[2].Y()-i_cond.GetGeometry()[0].Y())/std::abs((i_cond.GetGeometry()[2].Y()-i_cond.GetGeometry()[0].Y()));
                int vRow1 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[0].Y() + 1e-12*(direction_v));
                int vRow2 = FindKnotSpans1D(rSpansV, i_cond.GetGeometry()[2].Y() - 1e-12*(direction_v));

                const double direction_w = (i_cond.GetGeometry()[2].Z()-i_cond.GetGeometry()[0].Z())/std::abs((i_cond.GetGeometry()[2].Z()-i_cond.GetGeometry()[0].Z()));
                int wRow1 = FindKnotSpans1D(rSpansW, i_cond.GetGeometry()[0].Z() + 1e-12*(direction_w));
                int wRow2 = FindKnotSpans1D(rSpansW, i_cond.GetGeometry()[2].Z() - 1e-12*(direction_w));

                for (int v_row = std::min(vRow1, vRow2); v_row < std::max(vRow1,vRow2)+1; v_row++) {
                    for (int w_row = std::min(wRow1, wRow2); w_row < std::max(wRow1,wRow2)+1; w_row++) {
                        perpendicular_conditions_per_u_direction[v_row][w_row].push_back(i_cond.GetGeometry()[0].X());
                    }
                }
            }
        }

        // Sort the vertical conditions by the x coordinate
        for (IndexType vRow = 0; vRow < rSpansV.size() - 1; ++vRow) {
            for (IndexType wRow = 0; wRow < rSpansW.size() - 1; ++wRow) {
                std::sort(perpendicular_conditions_per_u_direction[vRow][wRow].begin(), perpendicular_conditions_per_u_direction[vRow][wRow].end());  
            }
        }   
        // KRATOS_WATCH(perpendicular_conditions_per_u_direction)
        
        // Loop over each row and creation of the integration points
        for (IndexType j = 0; j < rSpansV.size() - 1; j++) {
            for (IndexType k = 0; k < rSpansW.size() - 1; k++) {

                int starting_u_index = 0;
                int nextSwitchKnotSpan;
                if (perpendicular_conditions_per_u_direction[j][k].size() > 0) {
                    // nextSwitchKnotSpan = (verticalConditionsPerRow[j][0]+1e-13)/meshSizes_uv[0];
                    nextSwitchKnotSpan = FindKnotSpans1D(rSpansU, perpendicular_conditions_per_u_direction[j][k][0]+1e-13);

                }
                else {
                    nextSwitchKnotSpan = rSpansU.size() - 1;
                }
                int perpendicular_condition_index = 1;

                if (isOuterLoopDefined){
                    if (perpendicular_conditions_per_u_direction[j][k].empty()) {continue;}
                    starting_u_index = FindKnotSpans1D(rSpansU, perpendicular_conditions_per_u_direction[j][k][0]+1e-13);
                    nextSwitchKnotSpan = FindKnotSpans1D(rSpansU, perpendicular_conditions_per_u_direction[j][k][1]+1e-13);
                    perpendicular_condition_index = 2;
                }

                for (IndexType i = starting_u_index; i < rSpansU.size() - 1; i++) {

                    if (i == nextSwitchKnotSpan) {
                        if (perpendicular_condition_index < perpendicular_conditions_per_u_direction[j][k].size()) {
                            i = FindKnotSpans1D(rSpansU, perpendicular_conditions_per_u_direction[j][k][perpendicular_condition_index]+1e-13);
                        } else {
                            break;
                        }
                        perpendicular_condition_index++;
                        if (perpendicular_condition_index < perpendicular_conditions_per_u_direction[j][k].size()) {
                            nextSwitchKnotSpan = FindKnotSpans1D(rSpansU, perpendicular_conditions_per_u_direction[j][k][perpendicular_condition_index]+1e-13);
                        }
                        else {
                            nextSwitchKnotSpan = rSpansU.size() - 1;
                        }
                        perpendicular_condition_index++;
                    }            

                    const IndexType number_of_integration_points = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) 
                                                                    * rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1) 
                                                                    * rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(2);

                    IndexType initial_integration_size = rIntegrationPoints.size();


                    if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                        rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                    }
                    typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                    advance(integration_point_iterator, initial_integration_size);
                    
                    IntegrationPointUtilities::IntegrationPoints3D(
                        integration_point_iterator,
                        rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), 
                        rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1),
                        rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(2),
                        rSpansU[i], rSpansU[i + 1],
                        rSpansV[j], rSpansV[j + 1],
                        rSpansW[k], rSpansW[k + 1]);
                }
            } 
        }
    };

} // namespace Kratos.
