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
        mControlVarsUpdate = ZeroVector(mTotalNumControlVars);
        mSearchDirection = ZeroVector(mTotalNumControlVars);
        mProjection = ZeroVector(mTotalNumControlVars);
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

        int opt_itr = mrSettings["opt_itr"].GetInt();

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

        Vector mObjectiveGradients = ZeroVector(mTotalNumControlVars);
        double objective_value = 0.0;
        for(auto& objetive : mrSettings["objectives"]){
            objective_value = objetive["value"].GetDouble();
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
                L2_norm = std::sqrt(L2_norm);         

                int node_index=0;
                for(auto& node : objective_controlled_obj_model_part.Nodes()){
                    if (objective_controlled_object_size==3){
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(objective_control_gradient_name));
                        for(int i=0; i<objective_controlled_object_size; i++)
                            mObjectiveGradients[objective_controlled_object_start_index+objective_controlled_object_size*node_index+i] = nodal_gradient(i) * ((weight/mSumObjectiveWeights) * (1.0/L2_norm));
                    }
                    else if (objective_controlled_object_size==1) {
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(objective_control_gradient_name));
                        mObjectiveGradients[objective_controlled_object_start_index+node_index] = nodal_gradient * ((weight/mSumObjectiveWeights) * (1.0/L2_norm));
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
        for(auto& constraint : mrSettings["constraints"]){
            if(constraint["is_active"].GetBool()){
                active_constraints_violations[constraint_index] = std::pow((constraint["value"].GetDouble()-constraint["ref_value"].GetDouble())/constraint["ref_value"].GetDouble(),1.0) ;                
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
                    L2_norm = std::sqrt(L2_norm);
                    // now do the l2_norm scaling
                    for (int i=0;i<constraint_controlled_object_size*constraint_controlled_obj_model_part.Nodes().size();i++)
                        active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+i) *= (1.0/L2_norm);
                }
                constraint_index++;                
            }
        }

        for(int cons_i=0;cons_i<num_active_constraints;cons_i++)
            noalias(row(active_constraints_gradients,cons_i)) = row(active_constraints_gradients,cons_i)/norm_2(row(active_constraints_gradients,cons_i));

        if(num_active_constraints>0)
        {
            // compute feasible self.lin_solversearch direction
            Matrix NTN = prod(active_constraints_gradients,trans(active_constraints_gradients));
            Matrix I = IdentityMatrix(active_constraints_gradients.size1());
            Matrix NTN_inv(NTN.size1(), NTN.size2());

            mrSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

            Vector current_projection = - (mObjectiveGradients - prod(trans(active_constraints_gradients), Vector(prod(NTN_inv, Vector(prod(active_constraints_gradients, mObjectiveGradients))))));
            Vector current_correction = - prod(trans(active_constraints_gradients), Vector(prod(NTN_inv,active_constraints_violations)));
            // 

            double BB_scaling_projection = 1.0,BB_scaling_correction = 1.0;
            if(norm_2(mCorrection)>0.0 && opt_itr >1){
                BB_scaling_projection = abs(inner_prod(mControlVarsUpdate,mControlVarsUpdate)/inner_prod(mControlVarsUpdate,current_projection-mProjection));
                BB_scaling_correction = abs(inner_prod(mControlVarsUpdate,mControlVarsUpdate)/inner_prod(mControlVarsUpdate,current_correction-mCorrection));
            }

            current_projection *= BB_scaling_projection;
            current_correction *= BB_scaling_correction;             

            current_projection = 0.618 * current_projection + (1.0-0.618) * mProjection;
            current_correction = 0.618 * current_correction + (1.0-0.618) * mCorrection;

            mProjection = current_projection;
            mCorrection = current_correction;

            // current_projection *= BB_scaling_projection;
            // current_correction *= BB_scaling_correction; 

            // if(norm_2(current_correction)>norm_2(current_projection))
            //     current_correction *= norm_2(current_projection)/norm_2(current_correction);


            std::cout<<"norm_2(current_projection): "<<norm_2(current_projection)<<", norm_2(current_correction): "<<norm_2(current_correction)<<std::endl;
            // std::cout<<"BB_scaling_projection: "<<BB_scaling_projection<<", BB_scaling_correction: "<<BB_scaling_correction<<std::endl;
            // std::cout<<"mult 1 : "<<BB_scaling_projection * norm_2(current_projection)<<", mult 2: "<<BB_scaling_correction * norm_2(current_correction)<<std::endl;


            mSearchDirection = current_projection + current_correction;
            // mSearchDirection = current_projection + current_correction;
            
        }
        else{
            mSearchDirection = - mObjectiveGradients;
            mCorrection.clear();
            mProjection.clear();
        }
            
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

            mControlVarsUpdate = control_scaling_factor * mSearchDirection;

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

    Vector mControlVarsUpdate;
    Vector mSearchDirection;
    Vector mProjection;
    Vector mCorrection;
    LinearSolver<DenseSpace, DenseSpace>& mrSolver;
    int mTotalNumControlVars;
    double mSumObjectiveWeights;
    double alpha;


}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // GRADIENT_PROJECTION_H
