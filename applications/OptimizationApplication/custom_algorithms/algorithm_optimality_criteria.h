// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef OPTIMALITY_CRITERIA_H
#define OPTIMALITY_CRITERIA_H

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

class KRATOS_API(OPTIMIZATION_APPLICATION) AlgorithmOptimalityCriteria : public OptimizationAlgorithm
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(AlgorithmOptimalityCriteria);

    typedef Variable<double> VariableType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpace;

    AlgorithmOptimalityCriteria(std::string OptName, Model& rModel, LinearSolver<DenseSpace, DenseSpace>& rSolver, Parameters& rOptSettings)
    : OptimizationAlgorithm(OptName,"OPTIMALITY_CRITERIA",rModel,rOptSettings),mrSolver(rSolver)
    {
    }

    virtual ~AlgorithmOptimalityCriteria() {};

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
        mControlVars = ZeroVector(mTotalNumControlVars);
        mSearchDirection = ZeroVector(mTotalNumControlVars);

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
        mControlVars.clear();
        mSearchDirection.clear();

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
                int node_index=0;
                for(auto& node : objective_controlled_obj_model_part.Nodes()){
                    if (objective_controlled_object_size==3){
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(objective_control_gradient_name));
                        for(int i=0; i<objective_controlled_object_size; i++)
                            mObjectiveGradients[objective_controlled_object_start_index+objective_controlled_object_size*node_index+i] += nodal_gradient(i) * (weight/mSumObjectiveWeights);
                    }
                    else if (objective_controlled_object_size==1) {
                        const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(objective_control_gradient_name));
                        mObjectiveGradients[objective_controlled_object_start_index+node_index] += nodal_gradient * (weight/mSumObjectiveWeights);
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
                    for(auto& node : constraint_controlled_obj_model_part.Nodes()){
                        if (constraint_controlled_object_size==3){
                            const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(constraint_control_gradient_name));
                            for(int i=0; i<constraint_controlled_object_size; i++){
                                active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+constraint_controlled_object_size*node_index+i) = nodal_gradient(i);   
                            }
                        }
                        else if (constraint_controlled_object_size==1) {
                            const auto & nodal_gradient = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(constraint_control_gradient_name));
                            active_constraints_gradients(constraint_index,constraint_controlled_object_start_index+node_index) = nodal_gradient;
                        } 
                        node_index++;
                    }
                }
                constraint_index++;                
            }
        }

        // now fill control vars
        int index = 0;
        for(auto& control : mrSettings["controls"]){
            int control_size = control["size"].GetInt();
            auto control_var_name = control["variable_name"].GetString();
            for(auto& control_obj : control["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                for(auto& node : r_controlling_object.Nodes()){
                    if (control_size==3){
                        auto & nodal_control = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(control_var_name));
                        mControlVars[index] = nodal_control(0);
                        mControlVars[index+1] = nodal_control(1);
                        mControlVars[index+2] = nodal_control(2);
                        index += 3;
                    }
                    else if (control_size==1){
                        auto & nodal_control = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(control_var_name));
                        mControlVars[index] = nodal_control;
                        index += 1;
                    }
                }
            }
        }

        //now compute lagrange multipliers
        Vector lagrangs;
        if(num_active_constraints>0)
        {

            // compute feasible self.lin_solversearch direction
            Matrix NTN = prod(active_constraints_gradients,trans(active_constraints_gradients));
            Matrix I = IdentityMatrix(active_constraints_gradients.size1());
            Matrix NTN_inv(NTN.size1(), NTN.size2());

            mrSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

            lagrangs = Vector(prod(NTN_inv, Vector(prod(active_constraints_gradients, mObjectiveGradients))));
            std::cout<<"lagrangs : "<<lagrangs<<std::endl;
        }
        else
            KRATOS_ERROR << "AlgorithmOptimalityCriteria: there should be active constraints" << std::endl;


        // compute and set the updates
        index = 0;
        for(auto& control : mrSettings["controls"]){
            int control_size = control["size"].GetInt();
            auto control_update_name = control["update_name"].GetString();
            auto control_var_name = control["variable_name"].GetString();
            auto control_max_update = control["max_update"].GetDouble();
            for(auto& control_obj : control["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                for(auto& node : r_controlling_object.Nodes()){
                    if (control_size==3){
                        auto & nodal_update = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(control_update_name));
                        auto & nodal_coords = node.Coordinates();
                        auto & nodal_var = node.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(control_var_name));
                        Vector Di = ZeroVector(3);
                        for(int acons_i =0; acons_i<num_active_constraints;acons_i++){
                            Di(0) += std::abs(lagrangs(acons_i)*active_constraints_gradients(acons_i,index)/mObjectiveGradients(index));
                            Di(1) += std::abs(lagrangs(acons_i)*active_constraints_gradients(acons_i,index+1)/mObjectiveGradients(index+1));
                            Di(2) += std::abs(lagrangs(acons_i)*active_constraints_gradients(acons_i,index+2)/mObjectiveGradients(index+2));
                        }

                        nodal_update(0) = nodal_var(0) * (std::sqrt(Di(0))-1);
                        nodal_update(1) = nodal_var(1) * (std::sqrt(Di(1))-1);
                        nodal_update(2) = nodal_var(2) * (std::sqrt(Di(2))-1);
                        index += 3;
                    }
                    else if (control_size==1){
                        auto & nodal_update = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(control_update_name));
                        auto & nodal_var = node.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(control_var_name));
                        double Di = 0.0;
                        for(int acons_i =0; acons_i<num_active_constraints;acons_i++)
                            Di += std::abs(mObjectiveGradients(index)/(lagrangs(acons_i)*active_constraints_gradients(acons_i,index)));

                        double update = nodal_var * (std::sqrt(Di)-1);
                        if (update>control_max_update)
                            update = control_max_update;

                        if (update<-control_max_update)
                            update = -control_max_update;

                        nodal_update = update;

                        if ((update+nodal_var)<0)
                            nodal_update =-nodal_var;

                        if ((update+nodal_var)>1.0)
                            nodal_update =1.0-nodal_var;

                        // if (update+nodal_var)                        

                        // nodal_update = nodal_var * (std::sqrt(Di)-1);
                        // if ((nodal_update+nodal_var)<0.0)
                        //    nodal_update 
                        // nodal_var * (std::sqrt(Di)-1);
                        index += 1;
                    }
                }
            }
        }
    };



    Vector mObjectiveGradients;
    Vector mControlVarsUpdate;
    Vector mControlVars;
    Vector mSearchDirection;
    LinearSolver<DenseSpace, DenseSpace>& mrSolver;
    int mTotalNumControlVars;
    double mSumObjectiveWeights;


}; // Class OptimizationAlgorithm

}  // namespace Kratos.

#endif // OPTIMALITY_CRITERIA_H
