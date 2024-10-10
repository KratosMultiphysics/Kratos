// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Aron Noordam
//

// System includes

// External includes

// Project includes
#include "set_moving_load_process.h"

#include <utilities/function_parser_utility.h>
#include <utilities/mortar_utilities.h>

#include "utilities/interval_utility.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

SetMovingLoadProcess::SetMovingLoadProcess(ModelPart& rModelPart,
                                                            Parameters Settings)
                                                            : mrModelPart(rModelPart),
                                                            mParameters(Settings)
{
    Parameters default_parameters(R"(
        {
            "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "POINT_LOAD",
            "load"            : [0.0, 1.0, 0.0],
            "direction"       : [1,1,1],
            "velocity"        : 1,
            "origin"          : [0.0, 0.0, 0.0],
            "offset"          : 0.0
        }  )"
    );
    
    // Set default velocity as a string, if the input velocity is a string
    if (mParameters.Has("velocity")){
        if (mParameters["velocity"].IsString()){
            default_parameters["velocity"].SetString("1");
        }
    }

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);


    // check if load parameter has size 3
    KRATOS_ERROR_IF(mParameters["load"].size() != 3) <<
        "'load' has to have size 3!" << std::endl;

    // check if all elements in load parameter are either string or a number
    bool is_all_string = true;
    bool is_all_number = true;
    for (IndexType i = 0; i < mParameters["load"].size(); i++)
    {
        if (!mParameters["load"][i].IsString()){
            is_all_string = false;
        }
        if(!mParameters["load"][i].IsNumber()){
            is_all_number = false;
        }
    }

    KRATOS_ERROR_IF(!is_all_string && !is_all_number) << "'load' has to be a vector of numbers, or an array with strings" << std::endl;

}

std::vector<IndexType> SetMovingLoadProcess::FindNonRepeatingIndices(const std::vector<IndexType> IndicesVector)
{
    // Insert all array elements in hash
    // table
    std::unordered_map<IndexType, int> mp;
    for (IndexType i = 0; i < IndicesVector.size(); i++)
        mp[IndicesVector[i]]++;

    // Traverse through map only and
    std::vector<IndexType> non_repeating_indices;
    for (auto x : mp)
        if (x.second == 1)
            non_repeating_indices.push_back(x.first);

    return non_repeating_indices;
          
}


 bool SetMovingLoadProcess::IsSwapPoints(const double FirstCoord, const double SecondCoord, const int Direction)
{
    // swap points if points are sorted in opposite order compared to direction
    if ((FirstCoord < SecondCoord) && (Direction < 0)){
        return true;
    }
    if ((FirstCoord > SecondCoord) && (Direction > 0)){
        return true;
    }
    return false;
}


Condition& SetMovingLoadProcess::GetFirstConditionFromCoord(const double FirstCoord, const double SecondCoord, const int Direction, std::vector<Condition>& rEndConditions)
{
    // if center1 coord  < center2 coord and direction is positive
    if ((FirstCoord < SecondCoord) && (Direction > 0)){
        return rEndConditions[0];
    }
    // if center1 coord  > center2 coord and direction is negative
    if ((FirstCoord > SecondCoord) && (Direction < 0)){
        return rEndConditions[0];
    }
    return rEndConditions[1];
}

Condition& SetMovingLoadProcess::GetFirstCondition(const Point FirstPoint, const Point SecondPoint, const array_1d<int,3> Direction, std::vector<Condition>& rEndConditions)
{
    constexpr double tolerance = std::numeric_limits<double>::epsilon();

    // sort on x-coord, if x coords are equal, sort on y coord, if y coord is equal sort on z-coord
    if (std::abs(FirstPoint[0] - SecondPoint[0]) > tolerance){
        return GetFirstConditionFromCoord(FirstPoint[0], SecondPoint[0], Direction[0], rEndConditions);
    }
    if (std::abs(FirstPoint[1] - SecondPoint[1]) > tolerance){
        return GetFirstConditionFromCoord(FirstPoint[1], SecondPoint[1], Direction[1], rEndConditions);
    }
        return GetFirstConditionFromCoord(FirstPoint[2], SecondPoint[2], Direction[2], rEndConditions);
}



bool SetMovingLoadProcess::IsConditionReversed(const Condition& rCondition, const array_1d<int, 3> Direction)
{
    constexpr double tolerance = std::numeric_limits<double>::epsilon();

    auto& r_points = rCondition.GetGeometry().Points();
    if (std::abs(r_points[0].X0() - r_points[1].X0()) > tolerance){
        return IsSwapPoints(r_points[0].X0(), r_points[1].X0(), Direction[0]);
    }
    if (std::abs(r_points[0].Y0() - r_points[1].Y0()) > tolerance){
        return IsSwapPoints(r_points[0].Y0(), r_points[1].Y0(), Direction[1]);
    }
    return IsSwapPoints(r_points[0].Z0(), r_points[1].Z0(), Direction[2]);
}


std::vector<Condition> SetMovingLoadProcess::SortConditions(ModelPart::ConditionsContainerType& rUnsortedConditions, Condition& rFirstCondition)
{

    std::vector<Condition> unsorted_conditions_v(rUnsortedConditions.begin(), rUnsortedConditions.end());

    std::vector<Condition> sorted_conditions;
    std::vector<int> visited_indices;
    GeometricalObject::GeometryType& r_geom_first = rFirstCondition.GetGeometry();
    std::vector<IndexType> node_id_vector{ r_geom_first[0].Id(),r_geom_first[1].Id() };

    bool is_cond_reversed = mIsCondReversedVector[0];
    while (visited_indices.size() != unsorted_conditions_v.size()) {
        for (IndexType i = 0; i < unsorted_conditions_v.size(); i++) {
            Condition& r_cond = unsorted_conditions_v[i];
            GeometricalObject::GeometryType& r_geom = r_cond.GetGeometry();


            // check if current index is already added to sorted condition vector
            if (!std::count(visited_indices.begin(), visited_indices.end(), i)) {
                // check if geom has a shared node with previous geom
                if (std::find(node_id_vector.begin(), node_id_vector.end(), r_geom.Points()[0].Id()) != node_id_vector.end() || std::find(node_id_vector.begin(), node_id_vector.end(), r_geom.Points()[1].Id()) != node_id_vector.end()) {
                    if (sorted_conditions.size() == 0) {
                        // check if both nodes of geom are equal to nodes in start element, only do this to add the first element in the sorted conditions vector
                        if (std::find(node_id_vector.begin(), node_id_vector.end(), r_geom.Points()[0].Id()) != node_id_vector.end() && std::find(node_id_vector.begin(), node_id_vector.end(), r_geom.Points()[1].Id()) != node_id_vector.end()) {
                            node_id_vector = { r_geom[0].Id(),r_geom[1].Id() };
                            sorted_conditions.push_back(r_cond);
                            visited_indices.push_back(i);
                        }
                    } else {
                        // sort nodes in condition, such that new node is connected to previous condition
                        IndexType prev_id;
                        if (is_cond_reversed) {
                            prev_id = node_id_vector[0];
                        } else {
                            prev_id = node_id_vector[1];
                        }

                        if (prev_id != r_geom.Points()[0].Id()) {
                            is_cond_reversed = true;
                            mIsCondReversedVector.push_back(is_cond_reversed);
                        } else {
                            is_cond_reversed = false;
                            mIsCondReversedVector.push_back(is_cond_reversed);
                        }

                        // add condition to sorted conditions vector
                        node_id_vector = { r_geom[0].Id(),r_geom[1].Id() };
                        sorted_conditions.push_back(r_cond);
                        visited_indices.push_back(i);
                    }
                }
            }

        }
    }

    return sorted_conditions;
}

std::vector<Condition> SetMovingLoadProcess::FindEndConditions()
{
    std::vector<IndexType> node_id_vector;
    const array_1d<double, 3> origin_point = mParameters["origin"].GetVector();

    // get all end node ids ( not the middle nodes, in case of line3 conditions)
    // simultaneously check if origin point is on line
    bool condition_is_on_line = false;
    for (auto& r_cond : mrModelPart.Conditions()) {

        auto& r_geom = r_cond.GetGeometry();
        Point local_point;
        if (r_geom.IsInside(origin_point, local_point)){
            condition_is_on_line = true;
        }

        node_id_vector.push_back(r_geom[0].Id());
        node_id_vector.push_back(r_geom[1].Id());
    }

    // error if origin point is not on line
    KRATOS_ERROR_IF_NOT(condition_is_on_line) << "Origin point of moving load is not on line" << std::endl;

    // find non repeating node ids
    const std::vector<IndexType> non_repeating_node_ids = FindNonRepeatingIndices(node_id_vector);

    // error if model part does not have 1 end and 1 beginning
    KRATOS_ERROR_IF_NOT(non_repeating_node_ids.size() == 2) << "Moving load condition model part needs to be connected with a beginning and end" << std::endl;

    // find conditions at both ends of model part
    std::vector<Condition> end_conditions;
    for (Condition& r_cond : mrModelPart.Conditions()) {

        auto& r_geom = r_cond.GetGeometry();
        for (IndexType i = 0; i < r_geom.size(); i++){
            for (IndexType j = 0; j < non_repeating_node_ids.size(); j++){
                if (r_geom[i].Id() == non_repeating_node_ids[j]){
                    end_conditions.push_back(r_cond);
                }
            }
        }
    }
    return end_conditions;
}


void SetMovingLoadProcess::InitializeDistanceLoadInSortedVector()
{
    double global_distance = 0;
    // loop over sorted conditions
    for (IndexType i = 0; i < mSortedConditions.size(); ++i){
        auto& r_cond = mSortedConditions[i];
        auto& r_geom = r_cond.GetGeometry();
        const double element_length = r_geom.Length();

        Point local_point;

        // read origin point
        const array_1d<double, 3> origin_point = mParameters["origin"].GetVector();

        // if origin point is within the current condition, set the global distance of the load, else continue the loop
        if (r_geom.IsInside(origin_point, local_point)){
            const double local_to_global_distance = (local_point[0] + 1) / 2 * element_length;
            if (mIsCondReversedVector[i]){
                mCurrentDistance = global_distance + element_length - local_to_global_distance;
            } else {
                mCurrentDistance = global_distance + local_to_global_distance;
            }
            mCurrentDistance += mParameters["offset"].GetDouble();
            return;
        }

        // add element length of current condition to the global distance
        global_distance += element_length;
    }
}


void SetMovingLoadProcess::ExecuteInitialize()
{
    KRATOS_TRY
    if (!this->mrModelPart.GetProcessInfo()[IS_RESTARTED]){
        // clear load functions vector
        mLoadFunctions.clear();

        // check if load input is a function or numeric and add load to member variable
        if (mParameters["load"][0].IsString()) {
            mUseLoadFunction = true;
            for (IndexType i = 0; i < mParameters["load"].size(); ++i) {
                BasicGenericFunctionUtility load_function = BasicGenericFunctionUtility(mParameters["load"][i].GetString());
                mLoadFunctions.push_back(load_function);
            }
        }
        else {
            mUseLoadFunction = false;
        }

        // check if velocity input is a function or numeric and add velocity to member variable
        if (mParameters["velocity"].IsString()) {
            mUseVelocityFunction = true;
        }
        else {
            mUseVelocityFunction = false;
        }

        array_1d<int, 3> direction;

        for (IndexType i = 0; i < mParameters["direction"].size(); ++i) {
            direction[i] = mParameters["direction"][i].GetInt();
        }

        // get the two line condition elements at both sides of the model part
        std::vector<Condition> end_conditions = FindEndConditions();

        // find start condition 
        const Point center_1 = end_conditions[0].GetGeometry().Center();
        const Point center_2 = end_conditions[1].GetGeometry().Center();
        Condition& r_first_cond = GetFirstCondition(center_1, center_2, direction, end_conditions);

        // Initialise vector which indicates if nodes in condition are in direction of movement
        mIsCondReversedVector.clear();
        mIsCondReversedVector.push_back(IsConditionReversed(r_first_cond, direction));
        mSortedConditions = SortConditions(mrModelPart.Conditions(), r_first_cond);

        InitializeDistanceLoadInSortedVector();
    }
    KRATOS_CATCH("")
}


void SetMovingLoadProcess::ExecuteInitializeSolutionStep()
{

    array_1d<double, 3> load_vector;

    // retrieve load from load function if given
    if (mUseLoadFunction){
        // get current time
        const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);

        for (IndexType i =0; i< mLoadFunctions.size();++i){
            load_vector[i] = mLoadFunctions[i].CallFunction(0, 0, 0, current_time, 0, 0, 0);
        }
    } else {
        load_vector = mParameters["load"].GetVector();
    }
    
    double distance_cond = 0;

    // bool to check if load is already added, such that a load is not added twice if the load is exactly at a shared node.
    bool is_moving_load_added = false;

    // loop over sorted conditions vector
    for (IndexType i = 0; i < mSortedConditions.size(); ++i) {
        auto& r_cond = mSortedConditions[i];
        auto& r_geom = r_cond.GetGeometry();
        const double element_length = r_geom.Length();

        // if moving load is located at current condition element, apply moving load, else apply a zero load
        if (distance_cond + element_length >= mCurrentDistance && distance_cond <= mCurrentDistance && !is_moving_load_added){
            double local_distance;
            if (mIsCondReversedVector[i]){
                local_distance = distance_cond + element_length - mCurrentDistance;
            } else {
                local_distance = mCurrentDistance - distance_cond;
            }
            
            r_cond.SetValue(POINT_LOAD, load_vector);

            // distance is correct assuming nodes in condition are correctly sorted, the sorting is done while initializing this process
            r_cond.SetValue(MOVING_LOAD_LOCAL_DISTANCE, local_distance);
            is_moving_load_added = true;
        } else {
            r_cond.SetValue(POINT_LOAD, ZeroVector(3));
            r_cond.SetValue(MOVING_LOAD_LOCAL_DISTANCE, 0);
        }
        distance_cond += element_length;
    }    
}



void SetMovingLoadProcess::ExecuteFinalizeSolutionStep()
{
    double load_velocity;
    // retrieve load velocity from velocity function if given
    if (mUseVelocityFunction){
        // get current time
        const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);

        BasicGenericFunctionUtility velocity_function = BasicGenericFunctionUtility(mParameters["velocity"].GetString());

        // update velocity value
        load_velocity = velocity_function.CallFunction(0, 0, 0, current_time, 0, 0, 0);
    } else {
        load_velocity = mParameters["velocity"].GetDouble();
    }

    // move the load
    mCurrentDistance = mCurrentDistance + mrModelPart.GetProcessInfo().GetValue(DELTA_TIME) * load_velocity;
}

void SetMovingLoadProcess::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Process);
    rSerializer.save("SortedConditions", mSortedConditions);
    rSerializer.save("IsCondReversedVector", mIsCondReversedVector);
    rSerializer.save("UseLoadFunction", mUseLoadFunction);
    rSerializer.save("UseVelocityFunction", mUseVelocityFunction);
    rSerializer.save("CurrentDistance", mCurrentDistance);

}

void SetMovingLoadProcess::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Process);
    rSerializer.load("SortedConditions", mSortedConditions);
    rSerializer.load("IsCondReversedVector", mIsCondReversedVector);
    rSerializer.load("UseLoadFunction", mUseLoadFunction);
    rSerializer.load("UseVelocityFunction", mUseVelocityFunction);
    rSerializer.load("CurrentDistance", mCurrentDistance);
}

}  // namespace Kratos.
