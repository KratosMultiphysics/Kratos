//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

// System includes

// External includes

// Project includes
#include "set_moving_load_process.h"
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
            "variable_name"   : "MOVING_LOAD",
            "load"            : [0.0, 1.0, 0.0],
            "direction"       : [1,1,1],
            "velocity"        : 1
        }  )"
    );
    Parameters mParameters;
    //IntervalUtility interval_utility(mParameters);

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    KRATOS_ERROR_IF(mParameters["load"].GetVector().size() != 3) <<
        "'load' has to be a vector of doubles with size 3!" << std::endl;
}




std::vector<int> FindNonRepeatingIndices(std::vector<int> arr)
{
    // Insert all array elements in hash
    // table
    std::unordered_map<int, int> mp;
    for (int i = 0; i < arr.size(); i++)
        mp[arr[i]]++;

    // Traverse through map only and
    std::vector<int> non_repeating_indices;
    for (auto x : mp)
        if (x.second == 1)
            non_repeating_indices.push_back(x.first);

    return non_repeating_indices;
          
}

Condition& GetFirstConditionFromCoord(double first_coord, double second_coord, int direction, std::vector<Condition>& end_conditions)
{
    // if center1 x-coord  < center2 x-coord and direction is positive
    if ((first_coord < second_coord) && (direction > 0))
    {
        return end_conditions[0];
    }
    // if center1 x-coord  > center2 x-coord and direction is negative
    if ((first_coord > second_coord) && (direction < 0))
    {
        return end_conditions[0];
    }
    	return end_conditions[1];

}

Condition& GetFirstCondition(Point first_point, Point second_point, vector<int> direction, std::vector<Condition>& end_conditions)
{
    // sort on x-coord, if x coords are equal, sort on y coord, if y coord is equal sort on z-coord
    if (abs(first_point[0] - second_point[0]) > DBL_EPSILON)
    {
        return GetFirstConditionFromCoord(first_point[0], second_point[0], direction[0], end_conditions);
    }
    if (abs(first_point[1] - second_point[1]) > DBL_EPSILON)
    {
        return GetFirstConditionFromCoord(first_point[1], second_point[1], direction[1], end_conditions);
    }
        return GetFirstConditionFromCoord(first_point[2], second_point[2], direction[2], end_conditions);

}


std::vector<Condition> BubbleSortConditions(ModelPart::ConditionsContainerType& unsorted_conditions, Condition& first_condition)
{
    //ModelPart::ConditionsContainerType& sorted_conditions = unsorted_conditions;
    std::vector<Condition> unsorted_conditions_v(unsorted_conditions.begin(), unsorted_conditions.end());

    for (int i = 0; i < unsorted_conditions_v.size(); ++i)
    {
        int tmp = unsorted_conditions_v[i].Id();
        auto& r_geom = unsorted_conditions_v[i].GetGeometry();
        int tmp2 = r_geom.Points()[0].Id();
        int tmp3 = r_geom.Points()[1].Id();
        int tmp4 = 2;
        double a = r_geom[0].X();
        int b = 1 + 1;
    }


    std::vector<Condition> sorted_conditions;
    std::vector<int> visited_indices;
    GeometricalObject::GeometryType& r_geom_first = first_condition.GetGeometry();
    std::set<unsigned long long> node_id_vector{ r_geom_first[0].Id(),r_geom_first[1].Id() };

 
    while (visited_indices.size() != unsorted_conditions_v.size())
    {
        for (int i =0; i< unsorted_conditions_v.size(); i++)
        {
            Condition& r_cond = unsorted_conditions_v[i];
            GeometricalObject::GeometryType& r_geom = r_cond.GetGeometry();


            // check if current index is already added to sorted condition vector
            if (!std::count(visited_indices.begin(), visited_indices.end(), i))
            {
                // check if geom has a shared node with previous geom
                if (node_id_vector.find(r_geom.Points()[0].Id()) != node_id_vector.end() || node_id_vector.find(r_geom.Points()[1].Id()) != node_id_vector.end())
                {
                    if (sorted_conditions.size() == 0)
                    {
                        // check if both nodes of geom are equal to nodes in start element, only do this to add the first element in the sorted conditions vector
                        if (node_id_vector.find(r_geom[0].Id()) != node_id_vector.end() && node_id_vector.find(r_geom.Points()[1].Id()) != node_id_vector.end())
                        {
                            node_id_vector = { r_geom[0].Id(),r_geom[1].Id() };
                            sorted_conditions.push_back(r_cond);
                            visited_indices.push_back(i);
                        }
                    }
                    else
                    {
                        node_id_vector = { r_geom[0].Id(),r_geom[1].Id() };
                        sorted_conditions.push_back(r_cond);
                        visited_indices.push_back(i);
                    }
                }
            }
            
        }
    }

    // swap first 2 conditions if nessesary
    //r_geom = first_condition.GetGeometry();
    int tmp = first_condition.Id();
    for (int i = 0; i < sorted_conditions.size(); ++i)
    {
        int tmp2 = sorted_conditions[i].Id();
        double tmp12 = sorted_conditions[i].GetGeometry()[0].X();
        int tmp7 = 2;
    }
    

    return sorted_conditions;
}


void SetMovingLoadProcess::ExecuteInitialize()
{
    KRATOS_TRY
    mLoad = mParameters["load"].GetVector();
    const double current_time = mrModelPart.GetProcessInfo().GetValue(TIME);

    vector<int> direction = mParameters["direction"].GetVector();
    mLoadVelocity = mParameters["velocity"].GetDouble();
    mCurrentDistance = 0;

    std::vector<int> node_id_vector;
    std::vector<int> node_id_vector2;

    // get all end node ids ( not the middle nodes, in case of line3 conditions)
    for (auto& r_cond : mrModelPart.Conditions()) {

        auto geom = r_cond.GetGeometry();
        node_id_vector.push_back(geom[0].Id());
        node_id_vector.push_back(geom[1].Id());
    }

    // find non repeating node ids
    std::vector<int> non_repeating_node_ids = FindNonRepeatingIndices(node_id_vector);

    // error if model part does not have 1 end and 1 beginning
    KRATOS_ERROR_IF_NOT(non_repeating_node_ids.size() ==2) << "Moving load condition model part needs to be connected with a beginning and end" << std::endl;
    
    

    // find conditions at both ends of modelpart
    std::vector<Condition> end_conditions;
    for (Condition& r_cond : mrModelPart.Conditions()) {

        auto& r_geom = r_cond.GetGeometry();

        for (int i = 0; i < r_geom.size(); i++)
        {
            for (int j = 0; j < non_repeating_node_ids.size(); j++)
            {
                if (r_geom[i].Id() == non_repeating_node_ids[j])
                {
                    end_conditions.push_back(r_cond);
                }
            }    
        }
    }


    // find start condition 
    Point center_1 = end_conditions[0].GetGeometry().Center();
    Point center_2 = end_conditions[1].GetGeometry().Center();

    Condition& r_first_cond = GetFirstCondition(center_1, center_2, direction, end_conditions);

    mSortedConditions = BubbleSortConditions(mrModelPart.Conditions(), r_first_cond);

    auto tmp = r_first_cond.GetGeometry();

    int a = 1 + 1;
    //for (int i = 0; i < arr_size; i++) {
    //    node_id_vector[node_id_vector[i] % arr_size]
    //        = node_id_vector[node_id_vector[i] % arr_size] + arr_size;
    //}


    //IntervalUtility interval_utility(mParameters);
   /* if (interval_utility.IsInInterval(current_time)) {
        double total_area = 0.0;
        for (auto& r_cond : mrModelPart.Conditions()) {
            total_area += r_cond.GetGeometry().Area();
        }

        const double global_total_area = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(total_area);

        Vector force_by_area = mParameters["load"].GetVector() / global_total_area;

        for (auto& r_cond : mrModelPart.Conditions()) {
            const double area = r_cond.GetGeometry().Area();
            r_cond.SetValue(SURFACE_LOAD, force_by_area * area);
        }
    }*/
    KRATOS_CATCH("")
}

void SetMovingLoadProcess::ExecuteInitializeSolutionStep()
{
    double distance_cond = 0;


    // bool to check if load is already added, such that a load is not added twice if the load is exactly at a shared node.
    bool is_moving_load_added = false;

    for (auto& r_cond : mSortedConditions)
    {
        auto& r_geom = r_cond.GetGeometry();
        double element_length = r_geom.Length();

        // if moving load is located at current condition element, apply moving load, else apply a zero load
        if ((distance_cond + element_length >= mCurrentDistance) && (distance_cond <= mCurrentDistance) && !is_moving_load_added)
        {
            
            r_cond.SetValue(MOVING_LOAD, mLoad);

            // todo, currently distance is only correct when nodes are sorted in direction of distance
            r_cond.SetValue(MOVING_LOAD_LOCAL_DISTANCE, mCurrentDistance - distance_cond);
            is_moving_load_added = true;
        }
        else
        {
            r_cond.SetValue(MOVING_LOAD, ZeroVector(3));
        }
        distance_cond += element_length;
    }
 /*   const double current_time = mrModelPart.GetProcessInfo().GetValue(TIME);

    Vector direction = mParameters["direction"].GetVector();
    double velocity = mParameters["velocity"].GetDouble();



    IntervalUtility interval_utility(mParameters);
    if (interval_utility.IsInInterval(current_time)) {
        double total_area = 0.0;
        for (auto& r_cond : mrModelPart.Conditions()) {
            total_area += r_cond.GetGeometry().Area();
        }

        const double global_total_area = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(total_area);

        Vector force_by_area = mParameters["load"].GetVector() / global_total_area;

        for (auto& r_cond : mrModelPart.Conditions()) {
            const double area = r_cond.GetGeometry().Area();
            r_cond.SetValue(SURFACE_LOAD, force_by_area * area);
        }
    }*/
}

void SetMovingLoadProcess::ExecuteFinalizeSolutionStep()
{
    mCurrentDistance = mCurrentDistance + mrModelPart.GetProcessInfo().GetValue(DELTA_TIME) * mLoadVelocity;
}

}  // namespace Kratos.
