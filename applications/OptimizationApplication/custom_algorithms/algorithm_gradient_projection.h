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

#ifndef GRADIENT_PROJECTION_H
#define GRADIENT_PROJECTION_H

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
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "algorithm_base.h"

// ==============================================================================

namespace Kratos
{

class KRATOS_API(OPTIMIZATION_APPLICATION) AlgorithmGradientProjection : public OptimizationAlgorithm
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(AlgorithmGradientProjection);

    typedef Variable<double> VariableType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpace;

    AlgorithmGradientProjection(std::string OptName, Model& rModel, LinearSolver<DenseSpace, DenseSpace>& rSolver, Parameters& rOptSettings)
    : OptimizationAlgorithm(OptName,"gradient_projection",rModel,rOptSettings),mrSolver(rSolver)
    {
    }

    virtual ~AlgorithmGradientProjection() {};

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
        KRATOS_INFO("AlgorithmGradientProjection:Initialize: ") << " total number of design variables " <<mTotalNumControlVars<< std::endl;


        //resize the required vectors
        mSearchDirection = ZeroVector(mTotalNumControlVars);

        current_is_feasible = true;
        prev_is_feasible = true;
        sum_obj_improvement = 0.0;
        sum_obj_improvement_prev = 0.0;
        last_feasible_obj_val = 0.0;
        mSumObjectivesImprovements = 0.0;
        for(auto& objetive : mrSettings["objectives"]){
            double objective_improvement = objetive["objective_improvement"].GetDouble();
            mSumObjectivesImprovements +=  std::abs(objective_improvement);
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

        for(auto& constraint : mrSettings["constraints"]){
            constraint.AddEmptyArray("controlled_objects_start_index");
            constraint.AddEmptyArray("controlled_objects_size");
            for(long unsigned int constraint_control_i = 0; constraint_control_i<constraint["controls"].size(); constraint_control_i++){
                auto constraint_control_name = constraint["controls"][constraint_control_i].GetString();
                auto constraint_controlled_object_name = constraint["controlled_objects"][constraint_control_i].GetString();
                int constraint_control_object_start_index = 0;
                int constraint_control_object_size = 0;
                bool found = false;
                for(auto& control : mrSettings["controls"]){
                    auto control_name = control["name"].GetString();
                    int control_size = control["size"].GetInt();
                    for(auto& control_controlling_obj : control["controlling_objects"]){
                        auto control_controlling_obj_name =  control_controlling_obj.GetString();
                        ModelPart& control_controlling_obj_model_part = mrModel.GetModelPart(control_controlling_obj_name);
                        if (control_name==constraint_control_name && control_controlling_obj_name==constraint_controlled_object_name){
                            constraint_control_object_size = control_size;
                            found = true;
                            break;
                        }
                        constraint_control_object_start_index += control_size * control_controlling_obj_model_part.Nodes().size();
                    }
                    if(found)
                        break;
                }
                constraint["controlled_objects_start_index"].Append(constraint_control_object_start_index);
                constraint["controlled_objects_size"].Append(constraint_control_object_size);
            }
        }            
    };
    // --------------------------------------------------------------------------
    void CalculateSolutionStep() override {

        int opt_itr = mrSettings["opt_itr"].GetInt();
        int num_active_constraints = mrSettings["num_active_consts"].GetInt();

        Vector mObjectiveGradients = ZeroVector(mTotalNumControlVars);
        double objective_value = 0.0;
        double sum_objective_value = 0.0;
        sum_obj_improvement =0.0;
        for(auto& objetive : mrSettings["objectives"]){
            objective_value = objetive["value"].GetDouble();
            auto objective_prev_value = objetive["prev_itr_value"].GetDouble();
            double objective_weight = objetive["weight"].GetDouble();
            sum_objective_value += objective_value;
            double objective_improvement = objetive["objective_improvement"].GetDouble();

            if(opt_itr>1){
                double rel_change = 100 * (objective_value - objective_prev_value) / objective_prev_value;
                sum_obj_improvement += rel_change;
                double scale = 1.0;
                if(rel_change<0.0){
                    double ratio = objective_improvement/abs(rel_change);
                    if(ratio>1.2)
                        ratio = 1.2;
                    if(ratio<0.8)
                        ratio = 0.8;                    
                    scale = ratio;
                }
                else
                    scale = 1.2;
                
                objective_weight *= scale;
                objetive["weight"].SetDouble(objective_weight);
            }

        
            for(long unsigned int objective_control_i = 0; objective_control_i<objetive["controls"].size(); objective_control_i++){
                auto objective_control_i_name = objetive["controls"][objective_control_i].GetString();
                double objective_control_weight = 1.0;

                for(auto& control : mrSettings["controls"]){
                    auto control_name = control["name"].GetString();
                    if(control_name==objective_control_i_name)
                        objective_control_weight = control["objectives_weight"].GetDouble();
                }

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
                // L2_norm = 1.0;        

                int node_index=0;
                for(auto& node : objective_controlled_obj_model_part.Nodes()){
                    if (objective_controlled_object_size==3){
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(objective_control_gradient_name));
                        for(int i=0; i<objective_controlled_object_size; i++)
                            mObjectiveGradients[objective_controlled_object_start_index+objective_controlled_object_size*node_index+i] = nodal_gradient(i) * ((objective_weight*objective_control_weight) * (1.0/L2_norm));
                    }
                    else if (objective_controlled_object_size==1) {
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(objective_control_gradient_name));
                        mObjectiveGradients[objective_controlled_object_start_index+node_index] = nodal_gradient * ((objective_weight*objective_control_weight) * (1.0/L2_norm));
                    } 
                    node_index++;
                }             
            }
        }        
        mObjectiveGradients *= (1.0/norm_2(mObjectiveGradients));

        //now fill the active_constraints_gradients matrix
        int constraint_index = 0;
        Matrix active_constraints_gradients = ZeroMatrix(num_active_constraints,mTotalNumControlVars);
        Vector active_constraints_violations(num_active_constraints,0.0);
        current_is_feasible = true;
        for(auto& constraint : mrSettings["constraints"]){
            if(constraint["is_active"].GetBool()){
                auto const_type = constraint["type"].GetString();
                double current_violation = (constraint["value"].GetDouble()-constraint["ref_value"].GetDouble())/constraint["ref_value"].GetDouble();
                if(100 * std::abs(current_violation)>0.5)
                    current_is_feasible = false;
                double previous_violation = (constraint["prev_itr_value"].GetDouble()-constraint["ref_value"].GetDouble())/constraint["ref_value"].GetDouble();
                double relative_change = (constraint["value"].GetDouble()-constraint["prev_itr_value"].GetDouble())/constraint["prev_itr_value"].GetDouble();
                double weight = constraint["weight"].GetDouble();
                
                // analysis and set the weight
                if(opt_itr>1){
                    
                    if((const_type=="equality") || (const_type=="initial_value_equality")){
                        //first check the oscill
                        if((100 * std::abs(previous_violation)>1) && (100 * std::abs(current_violation)>1) && (((current_violation>0.0) && (previous_violation<0)) || ((current_violation<0.0) && (previous_violation>0))))
                            weight *= 0.95; 
                        else if((std::abs(current_violation)>std::abs(previous_violation)) && (100 * std::abs(previous_violation)>0.5))
                            weight *= 1.25;
                        else if((std::abs(current_violation)<std::abs(previous_violation)) && (100 * std::abs(current_violation)>0.5) && (100 * std::abs(relative_change)<1.0))
                            weight *= 1.25;
                        else if((100 * std::abs(previous_violation)<0.1) && (100 * std::abs(current_violation)<0.1))
                            weight *= 0.95;                                                                                
                    }
                    else{
                        bool prev_itr_is_active = constraint["prev_itr_is_active"].GetBool();
                        if(prev_itr_is_active){
                            if((std::abs(current_violation)>std::abs(previous_violation)) && (100 * std::abs(previous_violation)>1))
                                weight *= 1.25;                        
                            else if((std::abs(current_violation)<std::abs(previous_violation)) && (100 * std::abs(relative_change)<1.0))
                                weight *= 1.25;
                        }
                    }

                    if(weight<1.0)
                        weight = 1.0;                        
                } 
                constraint["weight"].SetDouble(weight);

                active_constraints_violations[constraint_index] = weight * (constraint["value"].GetDouble()-constraint["ref_value"].GetDouble()) / constraint["ref_value"].GetDouble();                
                for(long unsigned int constraint_control_i = 0; constraint_control_i<constraint["controls"].size(); constraint_control_i++){
                    auto constraint_control_i_name = constraint["controls"][constraint_control_i].GetString();
                    double constraint_control_weight = 1.0;

                    for(auto& control : mrSettings["controls"]){
                        auto control_name = control["name"].GetString();
                        if(control_name==constraint_control_i_name)
                            constraint_control_weight = control["constraints_weight"].GetDouble();
                    }

                    int constraint_controlled_object_start_index = constraint["controlled_objects_start_index"][constraint_control_i].GetInt();
                    int constraint_controlled_object_size = constraint["controlled_objects_size"][constraint_control_i].GetInt();
                    auto constraint_control_gradient_name = constraint["control_gradient_names"][constraint_control_i].GetString();
                    auto constraint_controlled_object_name = constraint["controlled_objects"][constraint_control_i].GetString();
                    ModelPart& constraint_controlled_obj_model_part = mrModel.GetModelPart(constraint_controlled_object_name);
                    int node_index=0;
                    double L2_norm =0;
                    for(auto& node : constraint_controlled_obj_model_part.Nodes()){
                        if (constraint_controlled_object_size==3){
                            const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(constraint_control_gradient_name));
                            for(int i=0; i<constraint_controlled_object_size; i++){
                                active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+constraint_controlled_object_size*node_index+i) = nodal_gradient(i); 
                                L2_norm += nodal_gradient(i) * nodal_gradient(i);     
                            }
                        }
                        else if (constraint_controlled_object_size==1) {
                            const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(constraint_control_gradient_name));
                            active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+node_index) = nodal_gradient;
                            L2_norm += nodal_gradient * nodal_gradient;
                        } 
                        node_index++;
                    }
                    // L2_norm = 1.0;
                    L2_norm = std::sqrt(L2_norm);
                    // now do the l2_norm scaling
                    for (long unsigned int i=0;i<constraint_controlled_object_size*constraint_controlled_obj_model_part.Nodes().size();i++)
                        active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+i) *= (constraint_control_weight/L2_norm);
                }
                constraint_index++;                
            }
        }

        for(int cons_i=0;cons_i<num_active_constraints;cons_i++)
            noalias(row(active_constraints_gradients,cons_i)) = row(active_constraints_gradients,cons_i)/norm_2(row(active_constraints_gradients,cons_i));        


        auto projection_step_size = mrSettings["projection_step_size"].GetDouble();
        auto correction_step_size = mrSettings["correction_step_size"].GetDouble();        

        // apply contraction if necessary
        if (opt_itr>2 && current_is_feasible && prev_is_feasible && (sum_objective_value>last_feasible_obj_val))
            mSumObjectivesImprovements *= 0.9;

        if (opt_itr>2 && current_is_feasible && prev_is_feasible && (sum_objective_value<last_feasible_obj_val))
            mSumObjectivesImprovements *= 1.01;            

        KRATOS_INFO("mSumObjectivesImprovements: ")<<mSumObjectivesImprovements<<std::endl;

        if(opt_itr>1){
        
            double scale = 1.0;
            double ratio = mSumObjectivesImprovements/abs(sum_obj_improvement);
            if(ratio>1.0 && sum_obj_improvement<0.0 && sum_obj_improvement_prev<0.0 && current_is_feasible && prev_is_feasible)
                scale = ratio; 

            if(ratio<1.0 && sum_obj_improvement<0.0 && sum_obj_improvement_prev<0.0 && current_is_feasible && prev_is_feasible)
                scale = ratio;

            if(scale>1.2)
                scale = 1.2; 
            if(scale<0.8) 
                scale = 0.8;                  

            if(sum_obj_improvement>0.0 && current_is_feasible && prev_is_feasible)
                scale = 0.8;

            projection_step_size *= scale;
            mrSettings["projection_step_size"].SetDouble(projection_step_size);

            if(current_is_feasible)
                last_feasible_obj_val = sum_objective_value;
        }

        sum_obj_improvement_prev = sum_obj_improvement;
        prev_is_feasible = current_is_feasible;         

         
        if(num_active_constraints>0)
        {
            // compute feasible self.lin_solversearch direction
            Matrix NTN = prod(active_constraints_gradients,trans(active_constraints_gradients));
            Matrix I = IdentityMatrix(active_constraints_gradients.size1());
            Matrix NTN_inv(NTN.size1(), NTN.size2());

            mrSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

            Vector projection = - (mObjectiveGradients - prod(trans(active_constraints_gradients), Vector(prod(NTN_inv, Vector(prod(active_constraints_gradients, mObjectiveGradients))))));
            Vector correction = - prod(trans(active_constraints_gradients), Vector(prod(NTN_inv,active_constraints_violations)));
            double current_sin_alpha = norm_2(projection);
            mrSettings["sin_alpha"].SetDouble(current_sin_alpha);
    
            mSearchDirection = (projection_step_size * projection / norm_2(projection)) + (correction_step_size * correction);
        }
        else
            mSearchDirection = - projection_step_size * mObjectiveGradients / norm_2(mObjectiveGradients);
            
        // set the updates
        int index = 0;
        for(auto& control : mrSettings["controls"]){
            int control_size = control["size"].GetInt();
            auto control_update_name = control["update_name"].GetString();

            for(auto& control_obj : control["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                for(auto& node : r_controlling_object.Nodes()){
                    if (control_size==3){
                        auto & nodal_update = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(control_update_name));
                        nodal_update(0) = mSearchDirection[index];
                        nodal_update(1) = mSearchDirection[index+1];
                        nodal_update(2) = mSearchDirection[index+2];
                        index += 3;
                    }
                    else if (control_size==1){
                        auto & nodal_update = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(control_update_name));
                        nodal_update = mSearchDirection[index];
                        index += 1;
                    }
                }
            }
        }
    };

    Vector mSearchDirection;
    LinearSolver<DenseSpace, DenseSpace>& mrSolver;
    int mTotalNumControlVars;
    double mSumObjectivesImprovements;
    double sum_obj_improvement;  
    double sum_obj_improvement_prev;
    double last_feasible_obj_val;
    bool current_is_feasible;
    bool prev_is_feasible;


}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // GRADIENT_PROJECTION_H
