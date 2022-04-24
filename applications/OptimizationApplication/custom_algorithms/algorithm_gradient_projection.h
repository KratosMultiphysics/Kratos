// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

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

        // get alpha
        alpha = mrSettings["alpha"].GetDouble();
        KRATOS_INFO("AlgorithmGradientProjection:Initialize: ") << " alpha value for correction is " <<alpha<< std::endl;

        //resize the required vectors
        mObjectiveGradients = ZeroVector(mTotalNumControlVars);
        mControlVarsUpdate = ZeroVector(mTotalNumControlVars);
        mSearchDirection = ZeroVector(mTotalNumControlVars);
        mCorrection = ZeroVector(mTotalNumControlVars);

        mSumObjectiveWeights = 0.0;
        for(auto& objetive : mrSettings["objectives"]){
            double weight = objetive["objective_weight"].GetDouble();
            mSumObjectiveWeights +=  std::abs(weight);
            objetive.AddEmptyArray("controlled_objects_start_index");
            objetive.AddEmptyArray("controlled_objects_size");
            for(int objective_control_i = 0; objective_control_i<objetive["controls"].size(); objective_control_i++){
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
            for(int constraint_control_i = 0; constraint_control_i<constraint["controls"].size(); constraint_control_i++){
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

        mObjectiveGradients.clear();
        mControlVarsUpdate.clear();
        mSearchDirection.clear();
        mCorrection.clear();

        // analyze the constraints
        int num_active_constraints = 0;
        for(auto& constraint : mrSettings["constraints"]){
            auto ref_value = constraint["ref_value"].GetDouble();
            auto value = constraint["value"].GetDouble();
            auto type = constraint["type"].GetString();
            Vector constraint_gradients = ZeroVector(mTotalNumControlVars);
            bool is_active = false;
            if(type == "equality")
                is_active = true;
            else if(type == "smaller_than" && value>ref_value)
                is_active = true;
            else if(type == "bigger_than" && value<ref_value)
                is_active = true;

            if(is_active)
                num_active_constraints++;
            
            if (!constraint.Has("is_active"))
                constraint.AddBool("is_active",is_active);
            else
                constraint["is_active"].SetBool(is_active);            
        }

       for(auto& objetive : mrSettings["objectives"]){
            double weight = objetive["objective_weight"].GetDouble();            
            for(int objective_control_i = 0; objective_control_i<objetive["controls"].size(); objective_control_i++){
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
                L2_norm = std::sqrt(L2_norm) * std::sqrt(objetive["controls"].size());

                std::cout<<"objective : "<<objetive["name"].GetString()<<", objective_controlled_object_name : "<<objective_controlled_object_name<<", L2_norm : "<<L2_norm<<std::endl;           

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

        //now fill the active_constraints_gradients matrix
        int constraint_index = 0;
        Matrix active_constraints_gradients = ZeroMatrix(num_active_constraints,mTotalNumControlVars);
        for(auto& constraint : mrSettings["constraints"]){
            if(constraint["is_active"].GetBool()){
                for(int constraint_control_i = 0; constraint_control_i<constraint["controls"].size(); constraint_control_i++){
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
                    L2_norm = std::sqrt(L2_norm) * std::sqrt(constraint["controls"].size());
                    std::cout<<"constraint : "<<constraint["name"].GetString()<<", objective_controlled_object_name : "<<constraint_controlled_object_name<<", L2_norm : "<<L2_norm<<std::endl; 
                    // now do the l2_norm scaling
                    for (int i=0;i<constraint_controlled_object_size*constraint_controlled_obj_model_part.Nodes().size();i++)
                        active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+i) *= (1.0/L2_norm);
                }
                constraint_index++;                
            }
        }


        if(num_active_constraints>0)
        {
            // compute feasible self.lin_solversearch direction
            Matrix NTN = prod(active_constraints_gradients,trans(active_constraints_gradients));
            Matrix I = IdentityMatrix(active_constraints_gradients.size1());
            Matrix NTN_inv(NTN.size1(), NTN.size2());

            mrSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

            mSearchDirection = - (mObjectiveGradients - prod(trans(active_constraints_gradients), Vector(prod(NTN_inv, Vector(prod(active_constraints_gradients, mObjectiveGradients))))));
            mSearchDirection /= norm_2(mSearchDirection);


            double sum_violations = 0.0;
            for(auto& constraint : mrSettings["constraints"]){
                if(constraint["is_active"].GetBool()){        
                    auto ref_value = constraint["ref_value"].GetDouble();
                    auto value = constraint["value"].GetDouble();
                    double violation_percentes = 100 * std::abs(value-ref_value)/std::abs(ref_value);
                    sum_violations += violation_percentes;
                }
            }

            for(auto& constraint : mrSettings["constraints"]){
                if(constraint["is_active"].GetBool()){        
                    Vector constraint_gradient = row(active_constraints_gradients,0);
                    auto ref_value = constraint["ref_value"].GetDouble();
                    auto value = constraint["value"].GetDouble();
                    double violation_percentes = 100 * std::abs(value-ref_value)/std::abs(ref_value);
                    double constraint_alpha = 0;
                    if(sum_violations>1)
                      constraint_alpha = (violation_percentes/sum_violations) * alpha;
                    else
                      constraint_alpha = (violation_percentes/sum_violations) * sum_violations * alpha;     
                    
                    std::cout<<" ****** Constraint "<<constraint["name"].GetString()<<", violation: "<<violation_percentes<<" %, constraint alpha: "<<constraint_alpha<<std::endl;
                    if (value>ref_value)
                        constraint_gradient *= -1;

                    mCorrection += constraint_alpha * constraint_gradient;
                }
            }
            if(sum_violations>1)
                mSearchDirection = (1.0-alpha) * mSearchDirection + mCorrection;
            else
                mSearchDirection = (1.0-alpha*sum_violations) * mSearchDirection + mCorrection;
        }
        else
            mSearchDirection = - mObjectiveGradients;

        

        mSearchDirection /= norm_2(mSearchDirection);


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
                for(auto& node : r_controlling_object.Nodes()){
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
    Vector mCorrection;
    LinearSolver<DenseSpace, DenseSpace>& mrSolver;
    int mTotalNumControlVars;
    double mSumObjectiveWeights;
    double alpha;


}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // GRADIENT_PROJECTION_H
