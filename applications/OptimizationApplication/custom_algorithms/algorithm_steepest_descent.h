//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#ifndef STEEPEST_DESCENT_H
#define STEEPEST_DESCENT_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "algorithm_base.h"

// ==============================================================================

namespace Kratos
{

class KRATOS_API(OPTIMIZATION_APPLICATION) AlgorithmSteepestDescent : public OptimizationAlgorithm
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(AlgorithmSteepestDescent);

    typedef Variable<double> VariableType;

    AlgorithmSteepestDescent(std::string OptName, Model& rModel, Parameters& rOptSettings)
    : OptimizationAlgorithm(OptName,"steepest_desceent",rModel,rOptSettings)
    {
    }

    virtual ~AlgorithmSteepestDescent() {};

    // --------------------------------------------------------------------------
    void Initialize() override {

        mTotalNumControlVars = 0;
        for(auto& control : mrSettings["controls"]){
            int control_num_vars = 0;
            int control_size = control["size"].GetInt();
            auto control_name = control["name"].GetString();
            for(auto& control_obj : control["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                control_num_vars += control_size * r_controlling_object.Nodes().size();
            }
            mTotalNumControlVars += control_num_vars;
        }

        //resize the required vectors
        mObjectiveGradients = ZeroVector(mTotalNumControlVars);
        mControlVarsUpdate = ZeroVector(mTotalNumControlVars);
        mSearchDirection = ZeroVector(mTotalNumControlVars);

        mSumObjectiveWeights = 0.0;
        for(auto& objetive : mrSettings["objectives"]){
            double weight = objetive["objective_weight"].GetDouble();
            mSumObjectiveWeights += std::abs(weight);
            objetive.AddEmptyArray("controlled_objects_start_index");
            objetive.AddEmptyArray("controlled_objects_size");
            for(long unsigned int objective_control_i = 0; objective_control_i<objetive["controls"].size(); objective_control_i++){
                auto objective_control_name = objetive["controls"][objective_control_i].GetString();
                auto objective_controlled_object_name = objetive["controlled_objects"][objective_control_i].GetString();
                int objective_control_object_start_index = 0;
                int objective_control_object_size = 0;
                bool found = false;
                for(auto& control : mrSettings["controls"]){
                    auto control_name = control["name"].GetString();
                    int control_size = control["size"].GetInt();
                    for(auto& control_controlling_obj : control["controlling_objects"]){
                        auto control_controlling_obj_name =  control_controlling_obj.GetString();
                        ModelPart& control_controlling_obj_model_part = mrModel.GetModelPart(control_controlling_obj_name);
                        if (control_name==objective_control_name && control_controlling_obj_name==objective_controlled_object_name){
                            objective_control_object_size = control_size;
                            found = true;
                            break;
                        }
                        objective_control_object_start_index += control_size * control_controlling_obj_model_part.Nodes().size();
                    }
                    if(found)
                        break;
                }
                objetive["controlled_objects_start_index"].Append(objective_control_object_start_index);
                objetive["controlled_objects_size"].Append(objective_control_object_size);
            }
        }             

    };
    // --------------------------------------------------------------------------
    void CalculateSolutionStep() override {

        mObjectiveGradients.clear();
        mControlVarsUpdate.clear();
        mSearchDirection.clear();

        for(auto& objetive : mrSettings["objectives"]){
            double weight = objetive["objective_weight"].GetDouble();
            
            for(long unsigned int objective_control_i = 0; objective_control_i<objetive["controls"].size(); objective_control_i++){
                int objective_controlled_object_start_index = objetive["controlled_objects_start_index"][objective_control_i].GetInt();
                int objective_controlled_object_size = objetive["controlled_objects_size"][objective_control_i].GetInt();
                auto objective_control_gradient_name = objetive["control_gradient_names"][objective_control_i].GetString();
                auto objective_controlled_object_name = objetive["controlled_objects"][objective_control_i].GetString();
                ModelPart& objective_controlled_obj_model_part = mrModel.GetModelPart(objective_controlled_object_name);
                double L2_norm =0;
                for(auto& node : objective_controlled_obj_model_part.Nodes()){
                    if (objective_controlled_object_size==3){
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(objective_control_gradient_name));
                        for(int i=0; i<objective_controlled_object_size; i++)                       
                            L2_norm += nodal_gradient(i) * nodal_gradient(i);
                    }
                    else if (objective_controlled_object_size==1) {
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(objective_control_gradient_name));
                        L2_norm += nodal_gradient * nodal_gradient;
                    }
                }
                L2_norm = std::sqrt(L2_norm);           

                int node_index=0;
                for(auto& node : objective_controlled_obj_model_part.Nodes()){
                    if (objective_controlled_object_size==3){
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(objective_control_gradient_name));
                        for(int i=0; i<objective_controlled_object_size; i++)
                            mObjectiveGradients[objective_controlled_object_start_index+objective_controlled_object_size*node_index+i] += nodal_gradient(i) * ((weight/mSumObjectiveWeights) * (1.0/L2_norm));
                    }
                    else if (objective_controlled_object_size==1) {
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(objective_control_gradient_name));
                        mObjectiveGradients[objective_controlled_object_start_index+node_index] += nodal_gradient * ((weight/mSumObjectiveWeights) * (1.0/L2_norm));
                    } 
                    node_index++;
                }             
            }
        }

        // compute search direction
        mSearchDirection = -mObjectiveGradients;

        // // compute control variables update
        // mControlVarsUpdate = 20 * mSearchDirection;

        // compute and set the updates
        int index = 0;
        for(auto& control : mrSettings["controls"]){
            int control_size = control["size"].GetInt();
            auto control_update_name = control["update_name"].GetString();
            auto control_max_update = control["max_update"].GetDouble();

            double control_max_abs_value = 0.0;
            int control_begin_index = index;

            for(auto& control_obj : control["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                long unsigned int node_counter = 0;
                while (node_counter<r_controlling_object.Nodes().size()){
                    if (control_size==3){
                        if (std::abs(mSearchDirection[control_begin_index])>control_max_abs_value)
                            control_max_abs_value = std::abs(mSearchDirection[control_begin_index]);
                        if (std::abs(mSearchDirection[control_begin_index+1])>control_max_abs_value)
                            control_max_abs_value = std::abs(mSearchDirection[control_begin_index+1]);
                        if (std::abs(mSearchDirection[control_begin_index+2])>control_max_abs_value)
                            control_max_abs_value = std::abs(mSearchDirection[control_begin_index+2]);
                        control_begin_index += control_size;
                    }
                    else if(control_size==1){
                        if (std::abs(mSearchDirection[control_begin_index])>control_max_abs_value)
                            control_max_abs_value = std::abs(mSearchDirection[control_begin_index]);
                        control_begin_index +=1;
                    }
                    node_counter++;
                }
            }

            double control_scaling_factor = control_max_update/control_max_abs_value;

            for(auto& control_obj : control["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                for(auto& node : r_controlling_object.Nodes()){
                    if (control_size==3){
                        auto & nodal_update = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(control_update_name));
                        nodal_update(0) = control_scaling_factor * mSearchDirection[index];
                        nodal_update(1) = control_scaling_factor * mSearchDirection[index+1];
                        nodal_update(2) = control_scaling_factor * mSearchDirection[index+2];
                        index += 3;
                    }
                    else if (control_size==1){
                        auto & nodal_update = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(control_update_name));
                        nodal_update = control_scaling_factor * mSearchDirection[index];
                        index += 1;
                    }
                }
            }
        }



    };



    Vector mObjectiveGradients;
    Vector mControlVarsUpdate;
    Vector mSearchDirection;
    int mTotalNumControlVars;
    double mSumObjectiveWeights;


}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // STEEPEST_DESCENT_H
