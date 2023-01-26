//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

// System includes

// External includes
#include "json/json_fwd.hpp"

// Project includes
#include "set_parameter_field_process.hpp"

#include <utilities/function_parser_utility.h>
#include <utilities/mortar_utilities.h>

#include "utilities/interval_utility.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

    SetParameterFieldProcess::SetParameterFieldProcess(ModelPart& rModelPart,
                                                            Parameters Settings)
                                                            : mrModelPart(rModelPart),
                                                            mParameters(Settings)
{
        std::cout << "test 1" << std::endl;

    // function type: python, cpp, input
    Parameters default_parameters(R"(
        {
            "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "CUSTOM",
            "func_type"       : "input",               
            "function"        : "100-2*y",
			"dataset"         : "dummy"
        }  )"
    );
    Parameters mParameters;

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

}


    void SetParameterFieldProcess::SetValueAtElement(Element::Pointer rElement, const Variable<double>& rVar, double Value)
    {

        // todo cleanup

        auto p_prop = rElement->pGetProperties();

        //Element::PropertiesType::Pointer new_prop(p_prop);
        Properties::Pointer new_prop = Kratos::make_shared<Properties>(*p_prop);
        //Properties::Pointer new_prop2;
        // Properties::Pointer new_prop(p_prop);
    	std::cout << "value before: " << p_prop->GetValue(rVar) << std::endl;
        std::cout << "value input: " << Value << std::endl;

        new_prop->SetValue(rVar, Value);

        std::cout << "value after: " << p_prop->GetValue(rVar) << std::endl;
        rElement->SetProperties(new_prop);
        std::cout << "value after2: " << rElement->GetProperties()[rVar] << std::endl;

        std::cout << "\n" << std::endl;
    }


void SetParameterFieldProcess::ExecuteInitialize()
{
    KRATOS_TRY

    if (!this->mrModelPart.GetProcessInfo()[IS_RESTARTED]){

        if (mParameters["func_type"].GetString().compare("input") == 0)
        {
            BasicGenericFunctionUtility parameter_function = BasicGenericFunctionUtility(mParameters["function"].GetString());

            //ModelPart::ElementsContainerType& elements = mrModelPart.Elements();

            const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);

            for (Element& r_element : mrModelPart.Elements()) {


                auto& r_geom = r_element.GetGeometry();
                const double val = parameter_function.CallFunction(r_geom.Center().X(), r_geom.Center().Y(), r_geom.Center().Z(), current_time, 0, 0, 0);
                const Variable<double>& r_var = KratosComponents< Variable<double> >::Get(mParameters["variable_name"].GetString());

                this->SetValueAtElement(&r_element, r_var, val);

            }
        }
        else if (mParameters["func_type"].GetString().compare("python") == 0)
        {

            
            auto dataset = mParameters["dataset"].GetString();
            const Variable<double>& r_var = KratosComponents< Variable<double> >::Get(mParameters["variable_name"].GetString());

            Parameters new_data = Parameters(dataset);

            std::vector<double> data_vector;
            for (Parameters& it : new_data)
            {
                data_vector.push_back(it.GetDouble());
            }

            IndexType i = 0;
            for (Element& r_element : mrModelPart.Elements())
            {
                this->SetValueAtElement(&r_element, r_var, data_vector[i]);
                i++;
            }
        }

    }
    KRATOS_CATCH("")
}

void SetParameterFieldProcess::ExecuteBeforeSolutionLoop()
    {
    int a = 1;
    // const Variable<double>& r_var = KratosComponents< Variable<double> >::Get(mParameters["variable_name"].GetString());
    //
    //     for (auto& r_element: mrModelPart.Elements())
    //     {
    //         auto& r_prop =  r_element.GetProperties();
    //
    //         double val = r_prop.GetValue(r_var);
    //
    //         auto& rProc = mrModelPart.GetProcessInfo();
    //         auto& r_geom = r_element.GetGeometry();
    //
    //         SizeType n_int_points = r_geom.IntegrationPointsNumber();
    //         std::vector<double> values_on_integration_points(n_int_points, val);
    //
    //         r_element.SetValuesOnIntegrationPoints(r_var, values_on_integration_points, rProc);
    //
    //     }

    //auto& r_prop = rElement->GetProperties();
    }

//
//void SetParameterFieldProcess::ExecuteInitializeSolutionStep()
//{
//
//    array_1d<double, 3> load_vector;
//
//    // retrieve load from load function if given
//    if (mUseTimeFunction){
//        // get current time
//        const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);
//
//        for (IndexType i =0; i< mLoadFunctions.size();++i){
//            load_vector[i] = mLoadFunctions[i].CallFunction(0, 0, 0, current_time, 0, 0, 0);
//        }
//    } else {
//        load_vector = mParameters["load"].GetVector();
//    }
//    
//    double distance_cond = 0;
//
//    // bool to check if load is already added, such that a load is not added twice if the load is exactly at a shared node.
//    bool is_moving_load_added = false;
//
//    // loop over sorted conditions vector
//    for (IndexType i = 0; i < mSortedConditions.size(); ++i) {
//        auto& r_cond = mSortedConditions[i];
//        auto& r_geom = r_cond.GetGeometry();
//        const double element_length = r_geom.Length();
//
//        // if moving load is located at current condition element, apply moving load, else apply a zero load
//        if ((distance_cond + element_length >= mCurrentDistance) && (distance_cond <= mCurrentDistance) && !is_moving_load_added){
//            double local_distance;
//            if (mIsCondReversedVector[i]){
//                local_distance = distance_cond + element_length - mCurrentDistance;
//            } else {
//                local_distance = mCurrentDistance - distance_cond;
//            }
//            
//            r_cond.SetValue(POINT_LOAD, load_vector);
//
//            // distance is correct assuming nodes in condition are correctly sorted, the sorting is done while initializing this process
//            r_cond.SetValue(MOVING_LOAD_LOCAL_DISTANCE, local_distance);
//            is_moving_load_added = true;
//        } else {
//            r_cond.SetValue(POINT_LOAD, ZeroVector(3));
//            r_cond.SetValue(MOVING_LOAD_LOCAL_DISTANCE, 0);
//        }
//        distance_cond += element_length;
//    }
//}
//
//
//
//void SetParameterFieldProcess::ExecuteFinalizeSolutionStep()
//{
//    double load_velocity;
//    // retrieve load velocity from velocity function if given
//    if (mUseVelocityFunction){
//        // get current time
//        const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);
//
//        BasicGenericFunctionUtility velocity_function = BasicGenericFunctionUtility(mParameters["velocity"].GetString());
//
//        // update velocity value
//        load_velocity = velocity_function.CallFunction(0, 0, 0, current_time, 0, 0, 0);
//    } else {
//        load_velocity = mParameters["velocity"].GetDouble();
//    }
//
//    // move the load
//    mCurrentDistance = mCurrentDistance + mrModelPart.GetProcessInfo().GetValue(DELTA_TIME) * load_velocity;
//}



}  // namespace Kratos.
