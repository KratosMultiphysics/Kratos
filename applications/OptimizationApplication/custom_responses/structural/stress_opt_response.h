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

#ifndef STRESS_OPT_RESPONSE_H
#define STRESS_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_responses/response.h"
#include "custom_elements/adjoint_small_displacement_element.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(OPTIMIZATION_APPLICATION) StressOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of StressOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(StressOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StressOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings, std::vector<LinearSolverType::Pointer>& rLinearSolvers )
        : Response(ResponseName,"stress",rModel, ResponseSettings){
            for(long unsigned int i=0;i<mrResponseSettings["control_types"].size();i++){
                auto control_type = mrResponseSettings["control_types"][i].GetString();
                if(control_type=="shape"){
                    std::string gradient_mode = mrResponseSettings["gradient_settings"]["gradient_mode"].GetString();
                    if (gradient_mode.compare("semi_analytic") == 0)
                    {
                        double delta = mrResponseSettings["gradient_settings"]["step_size"].GetDouble();
                        mDelta = delta;
                    }
                    else
                        KRATOS_ERROR << "StressOptResponse: Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;                    
                }
            }

            if(!mrResponseSettings.Has("yield_stress_type"))
                KRATOS_ERROR << "StressOptResponse: yield_stress_type should be provided in the stress response settings !! " << std::endl;

            mYieldStressType = mrResponseSettings["yield_stress_type"].GetString();
            if(!(mYieldStressType=="VON_MISES" || mYieldStressType=="PK2_COMP" || mYieldStressType=="PK2_TENS"))
                KRATOS_ERROR << "StressOptResponse: "<<mYieldStressType<<" is not supported and available yield_stress_types are VON_MISES, PK2_COMP and PK2_TENS !! "<< std::endl;


            if(!mrResponseSettings.Has("yield_stress_limit"))
                KRATOS_ERROR << "StressOptResponse: yield_stress_limit should be provided in the stress response settings !! " << std::endl; 
            mYieldStressLimit =  mrResponseSettings["yield_stress_limit"].GetDouble(); 

            if(!mrResponseSettings.Has("stress_penalty_factor"))
                KRATOS_ERROR << "StressOptResponse: stress_penalty_factor should be provided in the stress response settings !! " << std::endl; 
            mStressPenaltyFactor = mrResponseSettings["stress_penalty_factor"].GetDouble();

            if(!mrResponseSettings.Has("heaviside_penalty_factor"))
                KRATOS_ERROR << "StressOptResponse: heaviside_penalty_factor should be provided in the stress response settings !! " << std::endl; 
            mHeavisidePenaltyFactor = mrResponseSettings["heaviside_penalty_factor"].GetDouble();                            

            if(!mrResponseSettings.Has("heaviside_beta"))
                KRATOS_ERROR << "StressOptResponse: heaviside_beta should be provided in the stress response settings !! " << std::endl; 
            mBeta = mrResponseSettings["heaviside_beta"].GetDouble();   

            for(long unsigned int lin_i=0;lin_i<rLinearSolvers.size();lin_i++)
                rLinearSystemSolvers.push_back(rLinearSolvers[lin_i]);

        }

    /// Destructor.
    virtual ~StressOptResponse()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override {
        for(long unsigned int i=0;i<mrResponseSettings["evaluated_objects"].size();i++){
            auto eval_obj = mrResponseSettings["evaluated_objects"][i].GetString();
            ModelPart& eval_model_part = mrModel.GetModelPart(eval_obj);
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            auto control_type = mrResponseSettings["control_types"][i].GetString();

            KRATOS_ERROR_IF_NOT(eval_model_part.Elements().size()>0)
            <<"StressOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"StressOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;
           
        }        

        CreateAdjointComputeModelParts();
        CalculateNodeNeighbourCount();

        KRATOS_INFO("StressOptResponse: initialized by type: ")<<mYieldStressType<<", limit: "<<mYieldStressLimit<<", stress penalty: "<<mStressPenaltyFactor;
        KRATOS_INFO(", heaviside penalty: ")<<mHeavisidePenaltyFactor<<", beta: "<<mBeta<<std::endl;       

    };

    void CreateAdjointComputeModelParts(){

        for(long unsigned int i=0;i<mrResponseSettings["evaluated_objects"].size();i++){
            auto eval_obj = mrResponseSettings["evaluated_objects"][i].GetString();
            ModelPart& eval_model_part = mrModel.GetModelPart(eval_obj);

            // create adj model part
            ModelPart& root_model_part = eval_model_part.GetRootModelPart();
            std::string adj_model_part_name =  root_model_part.Name()+"_STRESS_ADJOINT";
            ModelPart* p_adj_model_part;
            p_adj_model_part = &(root_model_part.CreateSubModelPart(adj_model_part_name));
            mpADJModelParts.push_back(p_adj_model_part);

            // adding nodes
            for(auto& node : eval_model_part.Nodes())
                p_adj_model_part->AddNode(&node);

            // creating elements
            ModelPart::ElementsContainerType &rmesh_elements = p_adj_model_part->Elements();   

            for (int i = 0; i < (int)eval_model_part.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = eval_model_part.ElementsBegin() + i;
                Element::Pointer p_element = new AdjointSmallDisplacementElement(it->Id(), it->pGetGeometry(), eval_model_part.pGetElement(it->Id()));
                rmesh_elements.push_back(p_element);
            }             
        }

        // now add dofs and apply D BC
        for(long unsigned int model_i =0;model_i<mpADJModelParts.size();model_i++)
        {
            ModelPart* mpADJModePart = mpADJModelParts[model_i];
            for(auto& node_i : mpADJModePart->Nodes())
            {
                node_i.AddDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"));
                if(node_i.IsFixed(KratosComponents<Variable<double>>::Get("DISPLACEMENT_X"))){
                    node_i.Fix(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"));
                    auto& node_i_adj_dis = node_i.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"));
                    node_i_adj_dis = 0.0;
                }
                node_i.AddDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"));
                if(node_i.IsFixed(KratosComponents<Variable<double>>::Get("DISPLACEMENT_Y"))){
                    node_i.Fix(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"));
                    auto& node_i_adj_dis = node_i.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"));
                    node_i_adj_dis = 0.0;
                }
                node_i.AddDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"));
                if(node_i.IsFixed(KratosComponents<Variable<double>>::Get("DISPLACEMENT_Z"))){
                    node_i.Fix(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"));
                    auto& node_i_adj_dis = node_i.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"));
                    node_i_adj_dis = 0.0;
                }                                
            }
        }

        // create strategies
        for(long unsigned int model_i=0;model_i<mpADJModelParts.size();model_i++){

            StrategyType* mpStrategy = new StrategyType (*mpADJModelParts[model_i],rLinearSystemSolvers[model_i]);            
            mpStrategy->Initialize();
            mpStrategies.push_back(mpStrategy);
        }

    }

    void CalculateNodeNeighbourCount()
    {
        for(long unsigned int model_i =0;model_i<mpADJModelParts.size();model_i++)
        {
            ModelPart* mpADJModePart = mpADJModelParts[model_i];
            auto& r_nodes = mpADJModePart->Nodes();
            int mNumNodes = r_nodes.size();

            VariableUtils variable_utils;
            variable_utils.SetFlag(STRUCTURE,true,r_nodes);

            // Note: this should not be parallel, the operation is not threadsafe if the variable is uninitialized
            for (auto& r_node : r_nodes)
            {
                r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0);
            }

            mNumNodes = mpADJModePart->GetCommunicator().GetDataCommunicator().SumAll(mNumNodes);

            auto& r_elements = mpADJModePart->Elements();
            const int num_elements = r_elements.size();

            #pragma omp parallel for
            for (int i = 0; i < num_elements; i++)
            {
                auto i_elem = r_elements.begin() + i;
                auto& r_geom = i_elem->GetGeometry();
                for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
                {
                    auto& r_node = r_geom[i];
                    if (r_node.Is(STRUCTURE))
                    {
                        r_node.SetLock();
                        r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS) += 1;
                        r_node.UnSetLock();
                    }
                }
            }

            mpADJModePart->GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);

        }
    } 
    // --------------------------------------------------------------------------
    double CalculateElementStress(Element& elem_i, const ProcessInfo &rCurrentProcessInfo){        

        std::vector<double> gp_weights_vector;
        elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, rCurrentProcessInfo);        
        std::vector<double> stress_gp_vector;
 
        if(mYieldStressType=="PK2_TENS" || mYieldStressType=="PK2_COMP"){
            std::vector<Vector> PK2_stress_gp_vector;
            elem_i.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, PK2_stress_gp_vector, rCurrentProcessInfo);
            stress_gp_vector.resize(gp_weights_vector.size());
            for(IndexType i = 0; i < gp_weights_vector.size(); i++)
                stress_gp_vector[i] = PK2_stress_gp_vector[i][0] + PK2_stress_gp_vector[i][1] + PK2_stress_gp_vector[i][2]; 
        }            
        else if(mYieldStressType=="VON_MISES")
            elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), stress_gp_vector, rCurrentProcessInfo); 


        double elem_value = 0.0;
        for(IndexType i = 0; i < gp_weights_vector.size(); i++){
            double gp_stress = stress_gp_vector[i];
            double ratio = gp_stress/mYieldStressLimit;
            double pow_val = -2.0 * mBeta * (ratio-1);
            if(pow_val>700)
                pow_val = 700;
            if(pow_val<-700) 
                pow_val = -700;             
            elem_value += gp_weights_vector[i] * (1.0/(1+std::exp(pow_val))) * std::pow(ratio,mHeavisidePenaltyFactor);
        }

        return elem_value;          
    }
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        double intg_stress = 0.0;     
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const ProcessInfo &CurrentProcessInfo = r_eval_object.GetProcessInfo();
            // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
            for (auto& elem_i : r_eval_object.Elements())
            {
                const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
                if(element_is_active)
                    intg_stress += CalculateElementStress(elem_i,CurrentProcessInfo);
            }
        }
        return intg_stress;
    };    

    void ComputeAdjointRHS(ModelPart* mpADJModelPart){

        VariableUtils().SetHistoricalVariableToZero(ADJOINT_RHS, mpADJModelPart->Nodes());
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const ProcessInfo &rCurrentProcessInfo = r_eval_object.GetProcessInfo();
            for (auto& elem_i : r_eval_object.Elements())
            {

                std::vector<double> gp_weights_vector;
                elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, rCurrentProcessInfo);
                const unsigned int number_integration_points = gp_weights_vector.size();        
                std::vector<double> stress_gp_vector;
                std::vector<Matrix> stress_gp_tensor;
                std::vector<double> von_mises_sensitivity_prefactors(number_integration_points);               
        
                if(mYieldStressType=="PK2_TENS" || mYieldStressType=="PK2_COMP"){
                    std::vector<Vector> PK2_stress_gp_vector;
                    elem_i.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, PK2_stress_gp_vector, rCurrentProcessInfo);
                    stress_gp_vector.resize(gp_weights_vector.size());
                    for(IndexType i = 0; i < gp_weights_vector.size(); i++)
                        stress_gp_vector[i] = PK2_stress_gp_vector[i][0] + PK2_stress_gp_vector[i][1] + PK2_stress_gp_vector[i][2]; 
                }            
                else if(mYieldStressType=="VON_MISES"){
                    elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), stress_gp_vector, rCurrentProcessInfo); 
                    elem_i.CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, stress_gp_tensor, rCurrentProcessInfo);
                    for(IndexType k = 0; k < number_integration_points; ++k)
                    {
                        Matrix & PK2_k = stress_gp_tensor[k];
                        double radicant = 0.0;
                        radicant += PK2_k(0,0)*PK2_k(0,0) + PK2_k(1,1)*PK2_k(1,1) + PK2_k(2,2)*PK2_k(2,2);
                        radicant -= PK2_k(0,0)*PK2_k(1,1) + PK2_k(0,0)*PK2_k(2,2) + PK2_k(1,1)*PK2_k(2,2);
                        radicant += 3*PK2_k(0,1)*PK2_k(0,1) + 3*PK2_k(0,2)*PK2_k(0,2) + 3*PK2_k(1,2)*PK2_k(1,2);
                        von_mises_sensitivity_prefactors[k] = 0.5/std::sqrt(radicant);
                    }                    
                }
                    


                // Some working variables
                const SizeType num_nodes = elem_i.GetGeometry().PointsNumber();
                const SizeType dimension = elem_i.GetGeometry().WorkingSpaceDimension();
                const SizeType num_dofs_per_node = dimension;
                const SizeType num_dofs = num_nodes * num_dofs_per_node;
                

                // Build vector of variables containing the DOF-variables of the primal problem
                std::vector<const Variable<double>*> primal_solution_variable_list {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};

                // Build vector of variables containing the ADJOINT_RHS
                std::vector<const Variable<double>*> adj_rhs_variable_list {&ADJOINT_RHS_X, &ADJOINT_RHS_Y, &ADJOINT_RHS_Z};                

                // Store primal results and initialize deformation
                Vector initial_state_variables;
                initial_state_variables.resize(num_dofs, false);

                for (IndexType i = 0; i < num_nodes; ++i)
                {
                    const IndexType index = i * num_dofs_per_node;
                    for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
                    {
                        initial_state_variables[index + j] = elem_i.GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]);
                        elem_i.GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 0.0;
                    }
                }
                
                
                for (IndexType i = 0; i < num_nodes; ++i)
                {
                    for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
                    {
                        elem_i.GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 1.0;                        
                        double gp_sensitivity = 0.0;
                        for(IndexType k = 0; k < number_integration_points; ++k){
                            double gp_stress = stress_gp_vector[k];
                            double ratio = gp_stress/mYieldStressLimit;
                            double pow_val = -2.0 * mBeta * (ratio-1); 
                            if(pow_val>700)
                                pow_val = 700;
                            if(pow_val<-700) 
                                pow_val = -700;                            
                            
                            //derivative of heaviside
                            double heav_derv_mult = std::pow(ratio,mHeavisidePenaltyFactor);
                            heav_derv_mult *=  (1.0/(1+std::exp(pow_val))) * (1.0/(1+std::exp(pow_val)));
                            heav_derv_mult *= std::exp(pow_val);
                            heav_derv_mult *= 2 * mBeta;
                            heav_derv_mult /= mYieldStressLimit;

                            //derivative of penalty term
                            double penal_derv_mult = 0;
                            if(mHeavisidePenaltyFactor>0)
                                penal_derv_mult += (1.0/(1+std::exp(pow_val))) * (mHeavisidePenaltyFactor) * std::pow(ratio,mHeavisidePenaltyFactor-1) / mYieldStressLimit;

                            if(mYieldStressType=="PK2_TENS" || mYieldStressType=="PK2_COMP"){
                                std::vector<Vector> partial_stress_derivatives;
                                elem_i.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, partial_stress_derivatives, rCurrentProcessInfo);
                                gp_sensitivity += gp_weights_vector[k] * (heav_derv_mult + penal_derv_mult) * (partial_stress_derivatives[k][0] + partial_stress_derivatives[k][1] + partial_stress_derivatives[k][2]);
                                // std::cout<<"heav_derv_mult: "<<heav_derv_mult<<", penal_derv_mult: "<<penal_derv_mult<<", partial_stress_derivatives[k]: "<<partial_stress_derivatives[k]<<std::endl;
                            }                                
                            else if(mYieldStressType=="VON_MISES"){

                                std::vector<Matrix> partial_stress_derivatives;
                                elem_i.CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, partial_stress_derivatives, rCurrentProcessInfo);

                                double sensitivity_entry_k = 0.0;
                                Matrix & PK2_k = stress_gp_tensor[k];

                                sensitivity_entry_k += 2*PK2_k(0,0)*partial_stress_derivatives[k](0,0);
                                sensitivity_entry_k += 2*PK2_k(1,1)*partial_stress_derivatives[k](1,1);
                                sensitivity_entry_k += 2*PK2_k(2,2)*partial_stress_derivatives[k](2,2);

                                sensitivity_entry_k -= PK2_k(1,1)*partial_stress_derivatives[k](0,0);
                                sensitivity_entry_k -= PK2_k(0,0)*partial_stress_derivatives[k](1,1);

                                sensitivity_entry_k -= PK2_k(0,0)*partial_stress_derivatives[k](2,2);
                                sensitivity_entry_k -= PK2_k(2,2)*partial_stress_derivatives[k](0,0);

                                sensitivity_entry_k -= PK2_k(1,1)*partial_stress_derivatives[k](2,2);
                                sensitivity_entry_k -= PK2_k(2,2)*partial_stress_derivatives[k](1,1);

                                sensitivity_entry_k += 6*PK2_k(0,1)*partial_stress_derivatives[k](0,1);
                                sensitivity_entry_k += 6*PK2_k(0,2)*partial_stress_derivatives[k](0,2);
                                sensitivity_entry_k += 6*PK2_k(1,2)*partial_stress_derivatives[k](1,2);

                                sensitivity_entry_k *= von_mises_sensitivity_prefactors[k];
                                sensitivity_entry_k *= gp_weights_vector[k];
                                sensitivity_entry_k *= (heav_derv_mult + penal_derv_mult);
                                gp_sensitivity += sensitivity_entry_k;
                            }
                        }
                            
                        elem_i.GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 0.0;                            
                        elem_i.GetGeometry()[i].FastGetSolutionStepValue(*adj_rhs_variable_list[j]) += gp_sensitivity;
                    }
                }

                // Recall primal solution
                for (IndexType i = 0; i < num_nodes; ++i)
                {
                    const IndexType index = i * num_dofs_per_node;
                    for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
                        elem_i.GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_state_variables[index + j];
                }

            }
        }
    }

    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        // solve adjoints 
        for(long unsigned int model_i =0;model_i<mpADJModelParts.size();model_i++){
            ComputeAdjointRHS(mpADJModelParts[model_i]);
            mpStrategies[model_i]->Solve();
        }

        // compute adjoint-based sensitivities       
        for(long unsigned int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            const ProcessInfo &CurrentProcessInfo = controlled_model_part.GetProcessInfo();
            auto control_type = mrResponseSettings["control_types"][i].GetString();
            std::string grad_field_name;
            if(control_type=="shape"){
                grad_field_name = mrResponseSettings["gradient_settings"]["shape_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<array_1d<double,3>>>::Get(grad_field_name), controlled_model_part.Nodes());
            }                
            else if(control_type=="material"){
                grad_field_name = mrResponseSettings["gradient_settings"]["material_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<double>>::Get(grad_field_name), controlled_model_part.Nodes());                
            }

            //elems
			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="shape")
                        CalculateElementShapeGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="material")
                        CalculateElementMaterialGradients(elem_i,grad_field_name,CurrentProcessInfo);                                              
                }
            }

            //conds
			for (auto& cond_i : controlled_model_part.Conditions()){
				const bool cond_is_active = cond_i.IsDefined(ACTIVE) ? cond_i.Is(ACTIVE) : true;
				if(cond_is_active){
                    if(control_type=="shape")
                        CalculateConditionShapeGradients(cond_i,grad_field_name,CurrentProcessInfo);                                              
                }
            }

            //obj
			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="shape")
                        CalculateElementObjShapeGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="material")
                        CalculateElementObjMaterialGradients(elem_i,grad_field_name,CurrentProcessInfo);                                               
                }
            }



        }

		KRATOS_CATCH("");
 
    };

    void CalculateElementShapeGradients(Element& elem_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        Vector lambda;
        Vector RHS;

        // Get adjoint variables 
        GetElementAdjointVector(elem_i,lambda);

        // Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
        elem_i.CalculateRightHandSide(RHS, rCurrentProcessInfo);
        for (auto& node_i : elem_i.GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            Vector RHS_perturbed = Vector(RHS.size());
            Vector derived_RHS = Vector(RHS.size());
            
            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            elem_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = inner_prod(lambda, derived_RHS);


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            elem_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = inner_prod(lambda, derived_RHS);

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            elem_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = inner_prod(lambda, derived_RHS);

            // Assemble sensitivity to node
            noalias(node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name))) += gradient_contribution;
        }

    };

    void CalculateElementObjShapeGradients(Element& elem_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        double elem_before_pert_stress = CalculateElementStress(elem_i,rCurrentProcessInfo);

        for (auto& node_i : elem_i.GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            double pert_elem_stress = 0;
            
            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            pert_elem_stress = CalculateElementStress(elem_i,rCurrentProcessInfo);
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = (pert_elem_stress-elem_before_pert_stress) / mDelta;


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            pert_elem_stress = CalculateElementStress(elem_i,rCurrentProcessInfo);
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = (pert_elem_stress-elem_before_pert_stress) / mDelta;

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            pert_elem_stress = CalculateElementStress(elem_i,rCurrentProcessInfo);      
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = (pert_elem_stress-elem_before_pert_stress) / mDelta;

            // Assemble sensitivity to node
            noalias(node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name))) += gradient_contribution;
        }

    };

    void CalculateElementObjMaterialGradients(Element& elem_i, std::string material_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t number_of_nodes = r_this_geometry.size();     

        std::vector<double> gp_weights_vector;
        elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, rCurrentProcessInfo);        
        std::vector<double> stress_gp_vector;
 
        if(mYieldStressType=="PK2_TENS" || mYieldStressType=="PK2_COMP"){
            std::vector<Vector> PK2_stress_gp_vector;
            elem_i.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, PK2_stress_gp_vector, rCurrentProcessInfo);
            stress_gp_vector.resize(gp_weights_vector.size());
            for(IndexType i = 0; i < gp_weights_vector.size(); i++)
                stress_gp_vector[i] = PK2_stress_gp_vector[i][0] + PK2_stress_gp_vector[i][1] + PK2_stress_gp_vector[i][2]; 
        }            
        else if(mYieldStressType=="VON_MISES")
            elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), stress_gp_vector, rCurrentProcessInfo); 


        double elem_sens = 0.0;
        for(IndexType i = 0; i < gp_weights_vector.size(); i++){
            double gp_stress = stress_gp_vector[i];
            double ratio = gp_stress/mYieldStressLimit;
            double pow_val = -2.0 * mBeta * (ratio-1);
            if(pow_val>700)
                pow_val = 700;
            if(pow_val<-700) 
                pow_val = -700;             
            elem_sens += gp_weights_vector[i] * (1.0/(1+std::exp(pow_val))) * std::pow(ratio,mHeavisidePenaltyFactor);
        }


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pe_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PE_D_FD);
            r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(material_gradien_name)) += d_pe_d_fd * elem_sens / number_of_nodes;
        }


    };    

    void CalculateConditionShapeGradients(Condition& cond_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        Vector u;
        Vector lambda;
        Vector RHS;

        // Get state solution
        const auto& rConstCondRef = cond_i;
        rConstCondRef.GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        GetConditionAdjointVector(cond_i,lambda);

        cond_i.CalculateRightHandSide(RHS, rCurrentProcessInfo);
        if(u.size()==lambda.size() && RHS.size()==lambda.size()){

            for (auto& node_i : cond_i.GetGeometry())
            {
                array_3d gradient_contribution(3, 0.0);
                Vector RHS_perturbed = Vector(RHS.size());
                Vector derived_RHS = Vector(RHS.size());
                
                // x-direction
                node_i.GetInitialPosition()[0] += mDelta;
                node_i.Coordinates()[0] += mDelta;
                cond_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
                noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
                node_i.GetInitialPosition()[0] -= mDelta;
                node_i.Coordinates()[0] -= mDelta;
                gradient_contribution[0] = inner_prod(lambda, derived_RHS);


                // y-direction
                node_i.GetInitialPosition()[1] += mDelta;
                node_i.Coordinates()[1] += mDelta;
                cond_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
                noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
                node_i.GetInitialPosition()[1] -= mDelta;
                node_i.Coordinates()[1] -= mDelta;
                gradient_contribution[1] = inner_prod(lambda, derived_RHS);

                // z-direction
                node_i.GetInitialPosition()[2] += mDelta;
                node_i.Coordinates()[2] += mDelta;
                cond_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
                noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
                node_i.GetInitialPosition()[2] -= mDelta;
                node_i.Coordinates()[2] -= mDelta;
                gradient_contribution[2] = inner_prod(lambda, derived_RHS);

                // Assemble sensitivity to node
                noalias(node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name))) += gradient_contribution;
            }
        }

    };

    void CalculateElementMaterialGradients(Element& elem_i, std::string material_gradien_name, const ProcessInfo &rCurrentProcessInfo){


        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t number_of_nodes = r_this_geometry.size();

        Vector u;
        Vector lambda;
        
        // Get state solution
        const auto& rConstElemRef = elem_i;
        rConstElemRef.GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        GetElementAdjointVector(elem_i,lambda);


        Vector d_RHS_d_E;
        double e_curr = elem_i.GetProperties().GetValue(YOUNG_MODULUS);
        elem_i.GetProperties().SetValue(YOUNG_MODULUS,1.0);
        elem_i.CalculateRightHandSide(d_RHS_d_E,rCurrentProcessInfo);
        elem_i.GetProperties().SetValue(YOUNG_MODULUS,e_curr);


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pe_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PE_D_FD);
            r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(material_gradien_name)) += d_pe_d_fd * inner_prod(d_RHS_d_E,lambda) / number_of_nodes;
        }

    };    

    void GetElementAdjointVector(Element& elem_i, Vector &rAdjoints){

        const GeometryType &rgeom = elem_i.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = elem_i.GetGeometry().WorkingSpaceDimension();
        const unsigned int local_size = num_nodes * dimension;

        if (rAdjoints.size() != local_size)
            rAdjoints.resize(local_size, false);

        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        rAdjoints[index++] =
            rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"));
        rAdjoints[index++] =
            rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"));
        rAdjoints[index++] =
            rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"));
        }
    }

    void GetConditionAdjointVector(Condition& cond_i, Vector &rAdjoints){

        const GeometryType &rgeom = cond_i.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = cond_i.GetGeometry().WorkingSpaceDimension();
        const unsigned int local_size = num_nodes * dimension;

        if (rAdjoints.size() != local_size)
            rAdjoints.resize(local_size, false);

        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        rAdjoints[index++] =
            rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"));
        rAdjoints[index++] =
            rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"));
        rAdjoints[index++] =
            rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"));
        }
    }    


    // --------------------------------------------------------------------------
       
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "StressOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StressOptResponse";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<LinearSolverType::Pointer> rLinearSystemSolvers;
    std::vector<StrategyType*>mpStrategies;
    std::vector<ModelPart*> mpADJModelParts;

    // Initialized by class constructor

    
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{
    double mDelta;
    double mYieldStressLimit;
    double mBeta;
    double mStressPenaltyFactor;
    double mHeavisidePenaltyFactor;
    std::string mYieldStressType;
    

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      StressOptResponse& operator=(StressOptResponse const& rOther);

    /// Copy constructor.
//      StressOptResponse(StressOptResponse const& rOther);


    ///@}

}; // Class StressOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // STRESS_OPT_RESPONSE_H
