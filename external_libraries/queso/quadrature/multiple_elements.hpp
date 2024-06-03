// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MULITPLE_ELEMENTS_INCLUDE_HPP
#define MULITPLE_ELEMENTS_INCLUDE_HPP

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/element_container.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "queso/includes/parameters.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  QuadratureMultipleElements. Provides assembly opeartions for tensor-product quadrature rules of multiple non-trimmed elements.
 * @author Manuel Messmer
 * @brief  Provides assembly for 3D quadrature rules.
 * @details Available Quadrature rules:
 *          {GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2}
*/
template<typename TElementType>
class QuadratureMultipleElements {

public:
    ///@name Type Defintitions
    ///@{

    typedef TElementType ElementType;
    typedef ElementContainer<ElementType> ElementContainerType;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Assemble tensor product quadrature rules.
    /// @param rElements
    /// @param rNumberOfElements
    /// @param rIntegrationOrder
    /// @param Method Integration method
    static void AssembleIPs(ElementContainerType& rElements, const Vector3i& rNumberOfElements, const Vector3i& rIntegrationOrder, IntegrationMethodType Method) {
        // Loop over all 3 space dimensions
        // i = 0: x
        // i = 1: y
        // i = 2: z
        // Find neighbour relations!
        for( IndexType i = 0; i <= 2; ++i){
            bool local_end = false;
            bool found = false;
            IndexType current_id = 1;
            IndexType next_id = 0;
            std::vector<ElementType*> door_to_door_neighbours;
            door_to_door_neighbours.reserve(20);
            // Check if element with index=1 is part of rElements.
            const auto first_element = rElements.pGetElement(1, found);
            if( found && !first_element->IsTrimmed() )
                door_to_door_neighbours.push_back(first_element);

            if( rNumberOfElements[i] == 1){
                if(door_to_door_neighbours.size() > 0 ){
                        AssignNumberNeighbours(door_to_door_neighbours, i);
                }
                door_to_door_neighbours.clear();
            }
            std::size_t active_element_counter = 1;
            // Loop until all elements in rElements have beend visited/found
            while( active_element_counter < rElements.size() ){
                ElementType* neighbour = nullptr;

                if(i == 0)
                    neighbour = rElements.pGetNextElementInX(current_id, next_id, found, local_end);
                else if(i==1)
                    neighbour = rElements.pGetNextElementInY(current_id, next_id, found, local_end);
                else
                    neighbour = rElements.pGetNextElementInZ(current_id, next_id, found, local_end);

                if( found ){
                    active_element_counter++;
                    if(neighbour->IsTrimmed()){ // Is trimmed can only be evaluated if element is found (otherwise nullptr)
                        local_end = true;
                    }
                    else {
                        door_to_door_neighbours.push_back(neighbour);
                    }
                }

                // if( door_to_door_neighbours.size() >= 40 ){ // Maximum number of points is reached or element is trimmed
                //     local_end = true;
                // }

                if( local_end ){
                    if(door_to_door_neighbours.size() > 0 ){
                        AssignNumberNeighbours(door_to_door_neighbours, i);
                    }
                    door_to_door_neighbours.clear();
                }
                current_id = next_id;
            }
            if(door_to_door_neighbours.size() > 0 ){
                AssignNumberNeighbours(door_to_door_neighbours, i);
            }
            door_to_door_neighbours.clear();
        }

        // const auto element_it_begin_nei = rElements.begin();
        // const int number_neighbours = rElements.size();
        // for(int i = 0; i < number_neighbours; ++i){
        //     auto element_it = element_it_begin_nei + i;
        //     if( !(*element_it)->IsTrimmed()){
        //         auto local_lower_point = (*element_it)->GetLowerBoundParam();
        //         auto local_upper_point = (*element_it)->GetUpperBoundParam();
        //         double a = 0.5*(local_upper_point[0] + local_lower_point[0]);
        //         double b = 0.5*(local_upper_point[1] + local_lower_point[1]);
        //         double c = 0.5*(local_upper_point[2] + local_lower_point[2]);

        //         auto& points = (*element_it)->GetIntegrationPoints();
        //         points.push_back( IntegrationPoint(a, b, c, (*element_it)->NeighbourCoefficient() ) );
        //     }

        // }

        // Set all as not visited
        const auto element_it_begin = rElements.begin();
        for( IndexType i = 0; i < rElements.size(); ++i){
            auto element_it = element_it_begin + i;
            (*element_it)->SetVisited(false);
        }

        std::vector<ElementType*> box_neighbours{};
        box_neighbours.reserve(10000);
        typename std::vector<ElementType*>::iterator max_element_it{};
        ElementType* neighbour = nullptr;

        // Construct vector of sorted,unvisiteed elements.
        std::vector<ElementType*> sorted_univisited_elements{};
        sorted_univisited_elements.insert(sorted_univisited_elements.begin(), rElements.begin_to_ptr(), rElements.end_to_ptr());

        // Erase all trimmed elements.
        sorted_univisited_elements.erase(std::remove_if(sorted_univisited_elements.begin(), sorted_univisited_elements.end(), [](const auto& rValue) {
            return rValue->IsTrimmed(); }), sorted_univisited_elements.end());

        // Sort according to NeighbourCoefficient().
        std::sort( sorted_univisited_elements.begin(), sorted_univisited_elements.end(), [](auto &rLHs, auto &rRHS) -> bool
                    // If value is equal, sort according to ElementId(). This is done to keep same test results.
                    { if( std::abs((rLHs)->NeighbourCoefficient() - (rRHS)->NeighbourCoefficient()) < ZEROTOL){ return rLHs->GetId() <= rRHS->GetId();}
                    else { return (rLHs)->NeighbourCoefficient() > (rRHS)->NeighbourCoefficient(); }
                    });

        const std::array<int,6> direction_to_dimension = {0, 0, 1, 1, 2, 2};
        // int color_count = 1;
        bool stop = false;
        //int stop_count = 0;
        while( sorted_univisited_elements.size() > 0 && !stop ){

            max_element_it = sorted_univisited_elements.begin();
            (*max_element_it)->SetVisited(true);
            box_neighbours.push_back(*max_element_it);

            double max_neighbour_coefficient = 1.0;
            std::array<int,3> current_dimensions = {1, 1, 1};

            while( max_neighbour_coefficient > ZEROTOL && !stop ){
                std::array<double,6> neighbour_coeff = {0, 0, 0 ,0 ,0 ,0};
                std::array<std::vector<ElementType*>,6> tmp_neighbours{};

                // Check all four nighbours
                const auto element_it_begin = box_neighbours.begin();
                for( IndexType i = 0; i < box_neighbours.size(); ++i){
                    auto element_it = element_it_begin + i;
                    int current_id = (*element_it)->GetId();

                    for( int direction = 0; direction < 6; ++direction){
                        bool found = false;
                        if( !rElements.IsLast(current_id, direction) ){
                            neighbour = NextElement(rElements, current_id, found, direction);
                            if( found && !neighbour->IsVisited() && !neighbour->IsTrimmed() ){
                                neighbour_coeff[direction] += neighbour->NeighbourCoefficient();
                                tmp_neighbours[direction].push_back(neighbour);
                            }
                        }
                    }
                }
                bool valid_direction_found = false;
                while( !valid_direction_found && max_neighbour_coefficient > ZEROTOL && !stop){
                    int max_neighbour_coefficient_index = std::max_element(neighbour_coeff.begin(), neighbour_coeff.end()) - neighbour_coeff.begin();

                    max_neighbour_coefficient = neighbour_coeff[max_neighbour_coefficient_index];

                    int number_neighbours_in_max_direction = tmp_neighbours[max_neighbour_coefficient_index].size();

                    int dimension = direction_to_dimension[max_neighbour_coefficient_index];
                    int required_number_neighbours = 0;
                    if( dimension == 0)
                        required_number_neighbours = current_dimensions[1]*current_dimensions[2];
                    else if( dimension == 1)
                        required_number_neighbours = current_dimensions[0]*current_dimensions[2];
                    else if( dimension == 2)
                        required_number_neighbours = current_dimensions[0]*current_dimensions[1];

                    if( number_neighbours_in_max_direction ==  required_number_neighbours){
                        auto element_it_begin = tmp_neighbours[max_neighbour_coefficient_index].begin();
                        for(  IndexType i = 0; i < tmp_neighbours[max_neighbour_coefficient_index].size(); ++i){
                            auto element_it = element_it_begin + i;
                            (*element_it)->SetVisited(true);
                            box_neighbours.push_back(*element_it);
                            valid_direction_found = true;
                        }

                        current_dimensions[direction_to_dimension[max_neighbour_coefficient_index]]++;
                    }
                    else {
                        neighbour_coeff[max_neighbour_coefficient_index] = 0.0;
                    }
                }

            }

            StoreIntegrationPoints( box_neighbours, current_dimensions, rIntegrationOrder, Method);
            // const auto element_it_begin_nei = box_neighbours.begin();
            // const int number_neighbours = box_neighbours.size();
            // for(int i = 0; i < number_neighbours; ++i){
                // auto element_it = element_it_begin_nei + i;
                // if( !(*element_it)->IsTrimmed()){
                    // auto local_lower_point = (*element_it)->GetLowerBoundParam();
                    // auto local_upper_point = (*element_it)->GetUpperBoundParam();
                    // double a = 0.5*(local_upper_point[0] + local_lower_point[0]);
                    // double b = 0.5*(local_upper_point[1] + local_lower_point[1]);
                    // double c = 0.5*(local_upper_point[2] + local_lower_point[2]);

                    // auto& points = (*element_it)->GetIntegrationPoints();
                    // points.push_back( IntegrationPoint(a, b, c, color_count ) );
                // }

            // }
            // color_count++;

            box_neighbours.clear();
            // Erase all visiteed element to make search in next iteration faster.
            sorted_univisited_elements.erase(std::remove_if(sorted_univisited_elements.begin(), sorted_univisited_elements.end(), [](const auto& rValue) {
                    return rValue->IsVisited(); }), sorted_univisited_elements.end());
        }
    }

private:

    ///@todo Refactor this
    static void AssignNumberNeighbours(std::vector<ElementType*>& rElements, IndexType direction) {
        const auto element_it_begin = rElements.begin();
        const int number_neighbours = rElements.size();

        for(int i = 0; i < number_neighbours; ++i){
            auto element_it = element_it_begin + i;
            if( number_neighbours > 1)
                (*element_it)->SetNeighbourCoefficient(LinearFunction(i,number_neighbours), direction);
            else
                (*element_it)->SetNeighbourCoefficient(number_neighbours, direction);
        }
    }

    static ElementType* NextElement(ElementContainerType& rElements, std::size_t id, bool& found, int direction ) {
        bool dummy_local_end;
        std::size_t dummy_next_id;

        switch( direction )
        {
        case 0:
            return rElements.pGetNextElementInX(id, dummy_next_id, found, dummy_local_end);
        case 1:
            return rElements.pGetPreviousElementInX(id, dummy_next_id, found, dummy_local_end);
        case 2:
            return rElements.pGetNextElementInY(id, dummy_next_id, found, dummy_local_end);
        case 3:
            return rElements.pGetPreviousElementInY(id, dummy_next_id, found, dummy_local_end);
        case 4:
            return rElements.pGetNextElementInZ(id, dummy_next_id, found, dummy_local_end);
        case 5:
            return rElements.pGetPreviousElementInZ(id, dummy_next_id, found, dummy_local_end);
        default:
            QuESo_ERROR << "There are only 6 different directions.\n";
        }
    }

    static void StoreIntegrationPoints(std::vector<ElementType*>& rElements, std::array<int,3>& rNumberKnotspans,
            const Vector3i& rIntegrationOrder, IntegrationMethodType Method) {
        const auto element_it_begin = rElements.begin();
        // Find global extrem points (within box)
        PointType global_lower_point_param{1e10, 1e10, 1e10};
        PointType global_upper_point_param{-1e10, -1e10, -1e10};

        for( IndexType i = 0; i < rElements.size(); ++i){
            auto element_it = *(element_it_begin + i);
            const auto& lower_point = element_it->GetBoundsUVW().first;
            const auto& upper_point = element_it->GetBoundsUVW().second;
            if( lower_point[0] < global_lower_point_param[0] )
                global_lower_point_param[0] = lower_point[0];
            if( lower_point[1] < global_lower_point_param[1] )
                global_lower_point_param[1] = lower_point[1];
            if( lower_point[2] < global_lower_point_param[2] )
                global_lower_point_param[2] = lower_point[2];

            if( upper_point[0] > global_upper_point_param[0] )
                global_upper_point_param[0] = upper_point[0];
            if( upper_point[1] > global_upper_point_param[1] )
                global_upper_point_param[1] = upper_point[1];
            if( upper_point[2] > global_upper_point_param[2] )
                global_upper_point_param[2] = upper_point[2];
        }

        // Loop over all elements
        for( IndexType i = 0; i < rElements.size(); ++i){
            auto element_it = *(element_it_begin + i);

            // Local lower and upper points
            const auto lower_point_param = element_it->GetBoundsUVW().first;
            const auto upper_point_param = element_it->GetBoundsUVW().second;

            std::array<typename ElementType::IntegrationPoint1DVectorType, 3> tmp_integration_points{};

            for( int direction = 0; direction < 3; ++direction){
                const double distance_global = global_upper_point_param[direction] - global_lower_point_param[direction];
                const double length_global = std::abs(global_upper_point_param[direction] - global_lower_point_param[direction]);

                const auto p_integration_point_list =
                    IntegrationPointFactory1D::GetGGQ(rIntegrationOrder[direction], rNumberKnotspans[direction], Method);

                for( IndexType j = 0; j < p_integration_point_list->size(); ++j){
                    double position = global_lower_point_param[direction] + distance_global* (*p_integration_point_list)[j][0];
                    double weight = length_global *  (*p_integration_point_list)[j][1];
                    std::array<double, 2> tmp_point = {position, weight};
                    if( lower_point_param[direction]-EPS3 <= position && position < upper_point_param[direction]-EPS3){
                        if( lower_point_param[direction]+0.1*EPS3 > position) {
                            // Make sure point is clearly inside one element.
                            tmp_point[0] = lower_point_param[direction]+0.1*EPS3;
                        }
                        tmp_integration_points[direction].push_back(tmp_point);
                    }
                }

                const SizeType PointsInU = tmp_integration_points[0].size();
                const SizeType PointsInV = tmp_integration_points[1].size();
                const SizeType PointsInW = tmp_integration_points[2].size();

                for (SizeType u = 0; u < PointsInU; ++u) {
                    for (SizeType v = 0; v < PointsInV; ++v) {
                        for( SizeType w = 0; w < PointsInW; ++w) {
                            const double weight = tmp_integration_points[0][u][1]*tmp_integration_points[1][v][1]*tmp_integration_points[2][w][1];
                            element_it->GetIntegrationPoints().push_back(
                                                            IntegrationPoint( tmp_integration_points[0][u][0],
                                                                            tmp_integration_points[1][v][0],
                                                                            tmp_integration_points[2][w][0],
                                                                            weight ) );
                        }
                    }
                }
            }
        }
    }

    static bool AllElementsVisited(ElementContainerType& rElements) {
        const auto element_it_begin = rElements.begin();
        const int number_neighbours = rElements.size();
        for(int i = 0; i < number_neighbours; ++i){
            auto element_it = element_it_begin + i;
            if( !(*element_it)->IsTrimmed() && !(*element_it)->IsVisited() ){
                return false;
            }
        }

        return true;
    }

    static double LinearFunction(int x, int number_neighbours) {
        const double center = (double) (number_neighbours-1) / 2.0;
        const double delta = std::abs(center - x);

        double value = (1.0 - 0.9/center*delta)* (double)number_neighbours;

        QuESo_ERROR_IF(value < ZEROTOL) << "Value too low\n";

        return value;
    }
}; // End Class QuadratureMultipleElements

///@}

} // End namespace queso

#endif // MULITPLE_ELEMENTS_INCLUDE_HPP