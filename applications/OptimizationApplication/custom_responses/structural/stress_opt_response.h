// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef STRESS_OPT_RESPONSE_H
#define STRESS_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_responses/response.h"
#include "custom_elements/adjoint_solid_element.h"
#include "custom_elements/adjoint_shell_element.h"

// External includes
#include "custom_external_libraries/autodiff/forward/dual.hpp"
using namespace autodiff;

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

    struct StressComps
    {
        dual sxx;
        dual syy;
        dual szz;
        dual sxy;
        dual sxz;
        dual syz;
    };    
    
    /// Pointer definition of StressOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(StressOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StressOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings, LinearSolverType::Pointer pLinearSolver)
        : Response(ResponseName,"stress",rModel, ResponseSettings){
            for(int i=0;i<mrResponseSettings["control_types"].size();i++){
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
            if(!(mYieldStressType=="VON_MISES_STRESS" || mYieldStressType=="VON_MISES_STRESS_BOTTOM_SURFACE" || mYieldStressType=="VON_MISES_STRESS_TOP_SURFACE" \
                || mYieldStressType=="VON_MISES_STRESS_MIDDLE_SURFACE" || mYieldStressType=="MAX_PRINCIPAL_STRESS_TOP_SURFACE" || mYieldStressType=="MAX_PRINCIPAL_STRESS_MIDDLE_SURFACE" \ 
                || mYieldStressType=="MAX_PRINCIPAL_STRESS_BOTTOM_SURFACE" || mYieldStressType=="MAX_PRINCIPAL_STRESS" || mYieldStressType=="BENDING_VON_MISES"))
                KRATOS_ERROR << "StressOptResponse: "<<mYieldStressType<<" is not supported and available yield_stress_types are VON_MISES_STRESS, VON_MISES_STRESS_BOTTOM_SURFACE,\n"<< \
                "VON_MISES_STRESS_TOP_SURFACE, MAX_PRINCIPAL_STRESS_TOP_SURFACE, MAX_PRINCIPAL_STRESS_MIDDLE_SURFACE, MAX_PRINCIPAL_STRESS_BOTTOM_SURFACE, MAX_PRINCIPAL_STRESS, and VON_MISES_STRESS_MIDDLE_SURFACE !! "<< std::endl;

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

            mpLinearSystemSolver = pLinearSolver;

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

        BuiltinTimer timer;      

        //get the primal analysis model part
        auto primal_model_part_name = mrResponseSettings["analysis_model_part"].GetString();
        ModelPart& primal_model_part = mrModel.GetModelPart(primal_model_part_name);
        //first check whether it has elements
        KRATOS_ERROR_IF_NOT(primal_model_part.Elements().size()>0)
        <<"StressOptResponse::Initialize: analysis model part "<<primal_model_part_name<<" must have elements !"<<std::endl;
        //check if shell or solid problem
        if(primal_model_part.ElementsBegin()->GetGeometry().LocalSpaceDimension()>2)
            mIfShell = false;
        else      
            mIfShell = true;

        //evaluated objects must have element
        for(int i=0;i<mrResponseSettings["evaluated_objects"].size();i++){
            auto eval_obj = mrResponseSettings["evaluated_objects"][i].GetString();
            ModelPart& eval_model_part = mrModel.GetModelPart(eval_obj);
            KRATOS_ERROR_IF_NOT(eval_model_part.Elements().size()>0)
            <<"StressOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;
        } 

        //controlled objects must have element
        for(int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto control_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& control_model_part = mrModel.GetModelPart(control_obj);
            KRATOS_ERROR_IF_NOT(control_model_part.Elements().size()>0)
            <<"StressOptResponse::Initialize: controlled object "<<control_obj<<" must have elements !"<<std::endl;
        }                
    
        CreateAdjointComputeModelPart();
        CalculateNodeNeighbourCount();

        KRATOS_INFO("StressOptResponse: initialized by type: ")<<mYieldStressType<<", limit: "<<mYieldStressLimit<<", stress penalty: "<<mStressPenaltyFactor;
        KRATOS_INFO(", heaviside penalty: ")<<mHeavisidePenaltyFactor<<", beta: "<<mBeta<<std::endl;
        KRATOS_INFO("StressOptResponse::Initialize: ") << " Finished in " << timer.ElapsedSeconds() << " s." << std::endl;       

    };

    void CreateAdjointComputeModelPart(){

        //Create adjoint model part from primal 
        auto primal_model_part_name = mrResponseSettings["analysis_model_part"].GetString();
        ModelPart& primal_model_part = mrModel.GetModelPart(primal_model_part_name);
        ModelPart& root_model_part = primal_model_part.GetRootModelPart();
        std::string adj_model_part_name =  root_model_part.Name()+"_ADJOINT";
        if(root_model_part.HasSubModelPart(adj_model_part_name))
            mpADJModelPart = &(root_model_part.GetSubModelPart(adj_model_part_name));
        else{
            mpADJModelPart = &(root_model_part.CreateSubModelPart(adj_model_part_name));

            // adding nodes
            for(auto& node : primal_model_part.Nodes())
                mpADJModelPart->AddNode(&node);

            // adding elements
            ModelPart::ElementsContainerType &rmesh_elements = mpADJModelPart->Elements();   

            for (int i = 0; i < (int)primal_model_part.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = primal_model_part.ElementsBegin() + i;
                Element::Pointer p_element;
                if(mIfShell)
                    p_element = new AdjointShellElement(it->Id(), it->pGetGeometry(), primal_model_part.pGetElement(it->Id()));
                else
                    p_element = new AdjointSolidElement(it->Id(), it->pGetGeometry(), primal_model_part.pGetElement(it->Id()));

                rmesh_elements.push_back(p_element);
            }            
                
            // now add dofs and apply D BC
            for(auto& node_i : mpADJModelPart->Nodes())
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
                if(mIfShell){
                    node_i.AddDof(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_X"));
                    if(node_i.IsFixed(KratosComponents<Variable<double>>::Get("ROTATION_X"))){
                        node_i.Fix(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_X"));
                        auto& node_i_adj_rot = node_i.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_X"));
                        node_i_adj_rot = 0.0;
                    }
                    node_i.AddDof(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Y"));
                    if(node_i.IsFixed(KratosComponents<Variable<double>>::Get("ROTATION_Y"))){
                        node_i.Fix(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Y"));
                        auto& node_i_adj_dis = node_i.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Y"));
                        node_i_adj_dis = 0.0;
                    }
                    node_i.AddDof(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Z"));
                    if(node_i.IsFixed(KratosComponents<Variable<double>>::Get("ROTATION_Z"))){
                        node_i.Fix(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Z"));
                        auto& node_i_adj_dis = node_i.FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Z"));
                        node_i_adj_dis = 0.0;
                    }                    
                }                                
            }
        }

        // create strategy
        mpStrategy = new StrategyType (*mpADJModelPart,mpLinearSystemSolver);

    }

    void CalculateNodeNeighbourCount()
    {
        for(int i=0;i<mrResponseSettings["evaluated_objects"].size();i++){
            auto eval_obj = mrResponseSettings["evaluated_objects"][i].GetString();
            ModelPart& eval_model_part = mrModel.GetModelPart(eval_obj);

            auto& r_nodes = eval_model_part.Nodes();
            int mNumNodes = r_nodes.size();

            VariableUtils variable_utils;
            variable_utils.SetFlag(STRUCTURE,true,r_nodes);

            // Note: this should not be parallel, the operation is not threadsafe if the variable is uninitialized
            for (auto& r_node : r_nodes)
            {
                r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0);
            }

            mNumNodes = eval_model_part.GetCommunicator().GetDataCommunicator().SumAll(mNumNodes);

            auto& r_elements = eval_model_part.Elements();
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

            eval_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);            

        }
    } 
    // --------------------------------------------------------------------------
    double CalculateElementValue(Element& elem_i, const ProcessInfo &rCurrentProcessInfo){        

        std::vector<double> gp_weights_vector;
        std::vector<double> value_gp_vector;
        std::vector<StressComps> stress_gp_tensor;
        ComputeStressTensorAtGPs(elem_i,stress_gp_tensor,rCurrentProcessInfo);
        ComputeValuesAtGPs(elem_i,stress_gp_tensor,value_gp_vector,rCurrentProcessInfo);
        ComputeIntegrationWeightsAtGPs(elem_i,gp_weights_vector,rCurrentProcessInfo);
                     
        double elem_value = 0.0;
        for(IndexType i = 0; i < gp_weights_vector.size(); i++){
            double gp_stress = value_gp_vector[i];
            double ratio = gp_stress/mYieldStressLimit;
            double pow_val = -2.0 * mBeta * (ratio-1);
            if(pow_val>700)
                pow_val = 700;
            if(pow_val<-700) 
                pow_val = -700;

            if(mYieldStressType.find("BENDING") != std::string::npos)            
                elem_value += gp_weights_vector[i] * gp_stress;
            else
                elem_value += gp_weights_vector[i] * (1.0/(1+std::exp(pow_val))) * std::pow(ratio,mHeavisidePenaltyFactor);
        }

        return elem_value;          
    }
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        double intg_stress = 0.0; 
        double tot_vol = 0.0;    
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const ProcessInfo &CurrentProcessInfo = r_eval_object.GetProcessInfo();
            const std::size_t domain_size = r_eval_object.GetProcessInfo()[DOMAIN_SIZE];
            
            #pragma omp parallel for reduction(+:intg_stress,tot_vol)
            for (auto& elem_i : r_eval_object.Elements())
            {
                const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
                if(element_is_active){
                    intg_stress += CalculateElementValue(elem_i,CurrentProcessInfo);
                    if(mIfShell)
                        tot_vol += elem_i.GetGeometry().Area();
                    else
                        tot_vol += elem_i.GetGeometry().Volume();
                }
            }
        }
        return intg_stress/tot_vol;
    };

    void ComputeAggregationFunctionDerivative(const std::vector<double> & val_at_gp, const std::vector<double> & gp_int_w, std::vector<double> & der_at_gp){

        der_at_gp.resize(val_at_gp.size(),0.0);
        for(IndexType k = 0; k < val_at_gp.size(); ++k)
        {
            double ratio = val_at_gp[k]/mYieldStressLimit;
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

            if(mYieldStressType.find("BENDING") != std::string::npos)
                der_at_gp[k] = gp_int_w[k];
            else
                der_at_gp[k] = gp_int_w[k] * (heav_derv_mult + penal_derv_mult);
        }
    }

    void PrepareListsAndValues(Element& elem_i,std::vector<const Variable<double>*>& primal_list, std::vector<const Variable<double>*>& adj_rhs_list ){

        // Some working variables
        const SizeType num_nodes = elem_i.GetGeometry().PointsNumber();
        const SizeType dimension = elem_i.GetGeometry().WorkingSpaceDimension();
        SizeType num_dofs_per_node;
        if(mIfShell)
            num_dofs_per_node = 2 * dimension;
        else
            num_dofs_per_node = dimension;

        const SizeType num_dofs = num_nodes * num_dofs_per_node;  

        // Build vector of variables containing the DOF-variables of the primal problem
        primal_list.push_back(&DISPLACEMENT_X);
        primal_list.push_back(&DISPLACEMENT_Y);
        primal_list.push_back(&DISPLACEMENT_Z);
        if(mIfShell){
            primal_list.push_back(&ROTATION_X);
            primal_list.push_back(&ROTATION_Y);
            primal_list.push_back(&ROTATION_Z);
        }

        // Build vector of variables containing the ADJOINT_RHS
        adj_rhs_list.push_back(&ADJOINT_RHS_X);
        adj_rhs_list.push_back(&ADJOINT_RHS_Y);
        adj_rhs_list.push_back(&ADJOINT_RHS_Z);
        if(mIfShell){
            adj_rhs_list.push_back(&ADJOINT_RHS_ROT_X);
            adj_rhs_list.push_back(&ADJOINT_RHS_ROT_Y);
            adj_rhs_list.push_back(&ADJOINT_RHS_ROT_Z);
        }                

        // Store primal results and initialize deformation
        for (IndexType i = 0; i < num_nodes; ++i)
        {
            const IndexType index = i * num_dofs_per_node;
            for(IndexType j = 0; j < primal_list.size(); ++j)
                elem_i.GetGeometry()[i].FastGetSolutionStepValue(*primal_list[j]) = 0.0;
        }

    }

    void ComputeValuesAndDerivativesAtGPs(Element& elem_i, const std::vector<StressComps>& stress_gp_tensor, std::vector<double>& value_gp_vector, std::vector<Vector> & value_der_wrt_str_gp_vetor, const ProcessInfo &rCurrentProcessInfo){

        //compute values and dervs at GPs
        for(int gp_i=0;gp_i<stress_gp_tensor.size();gp_i++){
            if(mYieldStressType.find("VON_MISES") != std::string::npos){
                   StressComps k_stress_comps = stress_gp_tensor[gp_i];
                   dual value_at_gp = VonMisesAtGP(k_stress_comps); 
                   Vector der = ZeroVector(6); 
                   der[0] = derivative(VonMisesAtGP, wrt(k_stress_comps.sxx), at(k_stress_comps)); 
                   der[1] = derivative(VonMisesAtGP, wrt(k_stress_comps.syy), at(k_stress_comps));
                   der[2] = derivative(VonMisesAtGP, wrt(k_stress_comps.szz), at(k_stress_comps));
                   der[3] = derivative(VonMisesAtGP, wrt(k_stress_comps.sxy), at(k_stress_comps));
                   der[4] = derivative(VonMisesAtGP, wrt(k_stress_comps.sxz), at(k_stress_comps));
                   der[5] = derivative(VonMisesAtGP, wrt(k_stress_comps.syz), at(k_stress_comps));
                   value_gp_vector.push_back(double(value_at_gp));
                   value_der_wrt_str_gp_vetor.push_back(der);
            }
            else if(mYieldStressType.find("MAX_PRINCIPAL") != std::string::npos){
                   StressComps k_stress_comps = stress_gp_tensor[gp_i];
                   dual value_at_gp = MaxPrincipalAtGP(k_stress_comps); 
                   Vector der = ZeroVector(6); 
                   der[0] = derivative(MaxPrincipalAtGP, wrt(k_stress_comps.sxx), at(k_stress_comps)); 
                   der[1] = derivative(MaxPrincipalAtGP, wrt(k_stress_comps.syy), at(k_stress_comps));
                   der[2] = derivative(MaxPrincipalAtGP, wrt(k_stress_comps.szz), at(k_stress_comps));
                   der[3] = derivative(MaxPrincipalAtGP, wrt(k_stress_comps.sxy), at(k_stress_comps));
                   der[4] = derivative(MaxPrincipalAtGP, wrt(k_stress_comps.sxz), at(k_stress_comps));
                   der[5] = derivative(MaxPrincipalAtGP, wrt(k_stress_comps.syz), at(k_stress_comps));
                   value_gp_vector.push_back(double(value_at_gp));
                   value_der_wrt_str_gp_vetor.push_back(der);
            }                               
        }
    }

   void ComputeValuesAtGPs(Element& elem_i, const std::vector<StressComps>& stress_gp_tensor, std::vector<double>& value_gp_vector, const ProcessInfo &rCurrentProcessInfo){

        //compute values and dervs at GPs
        for(int gp_i=0;gp_i<stress_gp_tensor.size();gp_i++){
            if(mYieldStressType.find("VON_MISES") != std::string::npos){
                   StressComps k_stress_comps = stress_gp_tensor[gp_i];
                   dual value_at_gp = VonMisesAtGP(k_stress_comps); 
                   value_gp_vector.push_back(double(value_at_gp));
            } 
            else if(mYieldStressType.find("MAX_PRINCIPAL") != std::string::npos){
                   StressComps k_stress_comps = stress_gp_tensor[gp_i];
                   dual value_at_gp = MaxPrincipalAtGP(k_stress_comps);
                   value_gp_vector.push_back(double(value_at_gp));             
            }                              
        }
    }    

    void ComputeStressTensorAtGPs(Element& elem_i,std::vector<StressComps>& stress_gp_tensor, const ProcessInfo &rCurrentProcessInfo){

        std::vector<Matrix> tmp_stress_gp_tensor;
        //compute stress matrix at GPs
        if(mIfShell){
            if (mYieldStressType.find("MIDDLE_SURFACE") != std::string::npos)
                elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<Matrix>>::Get("SHELL_STRESS_MIDDLE_SURFACE"), tmp_stress_gp_tensor, rCurrentProcessInfo);
            if (mYieldStressType.find("TOP_SURFACE") != std::string::npos)
                elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<Matrix>>::Get("SHELL_STRESS_TOP_SURFACE"), tmp_stress_gp_tensor, rCurrentProcessInfo);
            if (mYieldStressType.find("BOTTOM_SURFACE") != std::string::npos)
                elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<Matrix>>::Get("SHELL_STRESS_BOTTOM_SURFACE"), tmp_stress_gp_tensor, rCurrentProcessInfo);          
            if (mYieldStressType.find("BENDING_VON_MISES") != std::string::npos){
                elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<Matrix>>::Get("SHELL_STRESS_TOP_SURFACE"), tmp_stress_gp_tensor, rCurrentProcessInfo);
                std::vector<Matrix> tmp_bot_stress_gp_tensor;
                elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<Matrix>>::Get("SHELL_STRESS_BOTTOM_SURFACE"), tmp_bot_stress_gp_tensor, rCurrentProcessInfo); 

                for(int k=0;k<tmp_stress_gp_tensor.size();k++){
                    noalias(tmp_stress_gp_tensor[k]) -= tmp_bot_stress_gp_tensor[k];
                    tmp_stress_gp_tensor[k] /= 2.0;
                }                
            }                                 
        }
        else
            elem_i.CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, tmp_stress_gp_tensor, rCurrentProcessInfo);


        // now fill 
        for(int k=0;k<tmp_stress_gp_tensor.size();k++){

            StressComps k_stress_comps;
            k_stress_comps.sxx = tmp_stress_gp_tensor[k](0,0);
            k_stress_comps.syy = tmp_stress_gp_tensor[k](1,1);
            k_stress_comps.szz = tmp_stress_gp_tensor[k](2,2);
            k_stress_comps.sxy = tmp_stress_gp_tensor[k](0,1);
            k_stress_comps.sxz = tmp_stress_gp_tensor[k](0,2);
            k_stress_comps.syz = tmp_stress_gp_tensor[k](1,2);
            stress_gp_tensor.push_back(k_stress_comps);
        }
    }

    void ComputeIntegrationWeightsAtGPs(Element& elem_i, std::vector<double>& gp_weights_vector, const ProcessInfo &rCurrentProcessInfo){

        if(mIfShell){
            std::vector<double> tmp_value_gp_vector;
            elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), tmp_value_gp_vector, rCurrentProcessInfo);
            gp_weights_vector.resize(tmp_value_gp_vector.size(),elem_i.GetGeometry().Area()/tmp_value_gp_vector.size());                
        }
        else
            elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, rCurrentProcessInfo);
    }        

    static dual VonMisesAtGP(const StressComps & stress_comps){
        dual vm = pow((stress_comps.sxx - stress_comps.syy),2);
        vm += pow((stress_comps.syy - stress_comps.szz),2);
        vm += pow((stress_comps.szz - stress_comps.sxx),2);
        vm *= 0.5;
        vm += 3 * pow((stress_comps.sxy),2);
        vm += 3 * pow((stress_comps.sxz),2);
        vm += 3 * pow((stress_comps.syz),2);
        vm = sqrt(vm);
        return vm;
    }    

    static dual MaxPrincipalAtGP(const StressComps & stress_comps){

        dual to_return = 0;
        dual I1, I2, I3;
        dual a, b, c;
        dual norm = pow(stress_comps.sxx,2) + pow(stress_comps.syy,2) + pow(stress_comps.szz,2);
        norm += pow(stress_comps.sxy,2) + pow(stress_comps.sxz,2) + pow(stress_comps.syz,2);
        norm = sqrt(norm);

        StressComps norm_stress_comps;
        norm_stress_comps.sxx = stress_comps.sxx;
        norm_stress_comps.syy = stress_comps.syy;
        norm_stress_comps.szz = stress_comps.szz;
        norm_stress_comps.sxy = stress_comps.sxy;
        norm_stress_comps.sxz = stress_comps.syz;
        norm_stress_comps.syz = stress_comps.sxz;
        
        static constexpr double tolerance = std::numeric_limits<double>::epsilon();
        if (norm < tolerance) norm = 1.0;

        norm_stress_comps.sxx /= norm;
        norm_stress_comps.syy /= norm;
        norm_stress_comps.szz /= norm;
        norm_stress_comps.sxy /= norm;
        norm_stress_comps.sxz /= norm;
        norm_stress_comps.syz /= norm;


        I1 = norm_stress_comps.sxx + norm_stress_comps.syy +norm_stress_comps.szz; 

        I2 = (norm_stress_comps.sxx + norm_stress_comps.szz) * norm_stress_comps.syy + norm_stress_comps.sxx * norm_stress_comps.szz +
        -norm_stress_comps.sxy * norm_stress_comps.sxy - norm_stress_comps.sxz * norm_stress_comps.sxz - norm_stress_comps.syz * norm_stress_comps.syz; 

        I3 = (norm_stress_comps.syy * norm_stress_comps.szz - norm_stress_comps.sxz * norm_stress_comps.sxz) * norm_stress_comps.sxx -
                norm_stress_comps.syy * norm_stress_comps.syz * norm_stress_comps.syz - norm_stress_comps.szz * norm_stress_comps.sxy * norm_stress_comps.sxy +
                2.0 * norm_stress_comps.sxy * norm_stress_comps.sxz * norm_stress_comps.syz;    

        dual II1 = pow(I1, 2);

        dual R = (2.0 * II1 * I1 - 9.0 * I2 * I1 + 27.0 * I3) / 54.0;
        dual Q = (3.0 * I2 - II1) / 9.0;             

        if (abs(Q) > tolerance) {
            dual cos_phi = R / (sqrt(-pow(Q, 3)));
            if (cos_phi >= 1.0)
                cos_phi = 1.0;
            else if (cos_phi <= -1.0)
                cos_phi = -1.0;
            dual phi = acos(cos_phi);
            dual phi_3 = phi / 3.0;

            dual aux1 = 2.0 * sqrt(-Q);
            dual aux2 = I1 / 3.0;
            dual deg_120 = 2.0 / 3.0 * 3.14159265359;             

            for (IndexType j = 0; j < 3; ++j){
                //if(norm * (aux2 + aux1 * cos(phi_3 + deg_120 * j))>to_return)
                    to_return += norm * (aux2 + aux1 * cos(phi_3 + deg_120 * j));
            } 
        } else {
                //if(stress_comps.sxx>to_return)
                    to_return += stress_comps.sxx;                
                //if(stress_comps.syy>to_return)
                    to_return += stress_comps.syy;
                //if(stress_comps.szz>to_return)
                    to_return += stress_comps.szz;
        }

        return to_return;

    }        

    void ComputeAdjointRHS(ModelPart* mpADJModelPart){        

        VariableUtils().SetHistoricalVariableToZero(ADJOINT_RHS, mpADJModelPart->Nodes());
        if(mIfShell)
            VariableUtils().SetHistoricalVariableToZero(ADJOINT_RHS_ROT, mpADJModelPart->Nodes());

        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const ProcessInfo &rCurrentProcessInfo = r_eval_object.GetProcessInfo();
            const std::size_t domain_size = r_eval_object.GetProcessInfo()[DOMAIN_SIZE];
            #pragma omp parallel for
            for (auto& elem_i : r_eval_object.Elements())
            {

                // make a copy of the element for following SA analysis
                auto& r_this_geometry = elem_i.GetGeometry();
                const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
                const std::size_t number_of_nodes = r_this_geometry.size();


                Element::NodesArrayType node_array;
                for (auto& node_i : elem_i.GetGeometry()){
                    Element::NodeType::Pointer node_p = node_i.Clone();
                    node_array.push_back(node_p);
                }

                Element::Pointer p_elem = elem_i.Create(elem_i.Id(),node_array, elem_i.pGetProperties());
                p_elem->SetData(elem_i.GetData());
                p_elem->Set(Flags(elem_i));
                p_elem->Initialize(rCurrentProcessInfo);                 

                std::vector<double> gp_weights_vector;
                std::vector<double> value_gp_vector;
                std::vector<Vector> value_der_wrt_str_gp_vetor;
                std::vector<double> sensitivity_prefactors;
                std::vector<double> aggregation_der_gp_vetor;
                std::vector<StressComps> stress_gp_tensor; 
                ComputeStressTensorAtGPs(*p_elem,stress_gp_tensor,rCurrentProcessInfo);
                ComputeValuesAndDerivativesAtGPs(*p_elem,stress_gp_tensor,value_gp_vector,value_der_wrt_str_gp_vetor,rCurrentProcessInfo);             
                ComputeIntegrationWeightsAtGPs(*p_elem,gp_weights_vector,rCurrentProcessInfo);
                ComputeAggregationFunctionDerivative(value_gp_vector,gp_weights_vector,aggregation_der_gp_vetor); 

                std::vector<const Variable<double>*> primal_solution_variable_list,adj_rhs_variable_list;
                PrepareListsAndValues(*p_elem,primal_solution_variable_list,adj_rhs_variable_list);                                          
                    

                // Some working variables
                const SizeType num_nodes = p_elem->GetGeometry().PointsNumber();
                const SizeType dimension = p_elem->GetGeometry().WorkingSpaceDimension();
                SizeType num_dofs_per_node;
                if(mIfShell)
                    num_dofs_per_node = 2 * dimension;
                else
                    num_dofs_per_node = dimension;

                const SizeType num_dofs = num_nodes * num_dofs_per_node;  
                                                                                                                      
                for (IndexType i = 0; i < num_nodes; ++i)
                {
                    const IndexType index = i * num_dofs_per_node;
                    for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
                    {
                        p_elem->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 1.0;                        
                        double gp_sensitivity = 0.0;
                        for(IndexType k = 0; k < aggregation_der_gp_vetor.size(); ++k){

                            std::vector<StressComps> partial_stress_derivatives;
                            ComputeStressTensorAtGPs(*p_elem,partial_stress_derivatives,rCurrentProcessInfo);

                            double agg_der = aggregation_der_gp_vetor[k];
                            double sensitivity_entry_k = 0;                            
                            sensitivity_entry_k += value_der_wrt_str_gp_vetor[k][0] * double(partial_stress_derivatives[k].sxx);
                            sensitivity_entry_k += value_der_wrt_str_gp_vetor[k][1] * double(partial_stress_derivatives[k].syy);
                            sensitivity_entry_k += value_der_wrt_str_gp_vetor[k][2] * double(partial_stress_derivatives[k].szz);
                            sensitivity_entry_k += value_der_wrt_str_gp_vetor[k][3] * double(partial_stress_derivatives[k].sxy);
                            sensitivity_entry_k += value_der_wrt_str_gp_vetor[k][4] * double(partial_stress_derivatives[k].sxz);
                            sensitivity_entry_k += value_der_wrt_str_gp_vetor[k][5] * double(partial_stress_derivatives[k].syz);
                            sensitivity_entry_k *= agg_der;
                            gp_sensitivity += sensitivity_entry_k;
                        }
                            
                        p_elem->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = 0.0;
                        #pragma omp atomic                            
                        elem_i.GetGeometry()[i].FastGetSolutionStepValue(*adj_rhs_variable_list[j]) += gp_sensitivity;
                    }
                }

            }
        }        
    }

    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        BuiltinTimer timer;

        // solve adjoints 
        ComputeAdjointRHS(mpADJModelPart);
        mpStrategy->Solve();

        KRATOS_INFO("StressOptResponse::CalculateGradient: ") << " Adjoint problem finished in " << timer.ElapsedSeconds() << " s." << std::endl;

        // compute adjoint-based sensitivities       
        for(int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
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
            else if(control_type=="thickness"){
                grad_field_name = mrResponseSettings["gradient_settings"]["thickness_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<double>>::Get(grad_field_name), controlled_model_part.Nodes());                
            }            


            //elems
            #pragma omp parallel for
			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="shape")
                        CalculateElementShapeGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="material")
                        CalculateElementMaterialGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="thickness")
                        CalculateElementThicknessGradients(elem_i,grad_field_name,CurrentProcessInfo);                                                                      
                }
            }

            //conds
			// for (auto& cond_i : controlled_model_part.Conditions()){
			// 	const bool cond_is_active = cond_i.IsDefined(ACTIVE) ? cond_i.Is(ACTIVE) : true;
			// 	if(cond_is_active){
            //         // if(control_type=="shape")
            //         //     CalculateConditionShapeGradients(cond_i,grad_field_name,CurrentProcessInfo);                                              
            //     }
            // }

            //obj
            #pragma omp parallel for
			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="shape")
                        CalculateElementObjShapeGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="material")
                        CalculateElementObjMaterialGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="thickness")
                        CalculateElementObjThicknessGradients(elem_i,grad_field_name,CurrentProcessInfo);                                                                       
                }
            }
        }

		KRATOS_CATCH("");
 
    };

    void CalculateElementShapeGradients(Element& elem_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();


        Element::NodesArrayType node_array;
        for (auto& node_i : elem_i.GetGeometry()){
            Element::NodeType::Pointer node_p = node_i.Clone();
            node_array.push_back(node_p);
        }

        Element::Pointer p_elem = elem_i.Create(elem_i.Id(),node_array, elem_i.pGetProperties());
        p_elem->SetData(elem_i.GetData());
        p_elem->Set(Flags(elem_i));
        p_elem->Initialize(rCurrentProcessInfo); 

        Vector lambda;
        Vector RHS;

        // Get adjoint variables 
        GetElementAdjointVector(*p_elem,lambda);

        // Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
        p_elem->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        int node_iter = 0;
        for (auto& node_i : p_elem->GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            Vector RHS_perturbed = Vector(RHS.size());
            Vector derived_RHS = Vector(RHS.size());
            
            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            p_elem->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = inner_prod(lambda, derived_RHS);


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            p_elem->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = inner_prod(lambda, derived_RHS);

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            p_elem->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = inner_prod(lambda, derived_RHS);

            // Assemble sensitivity to node
            array_3d& r_nodal_variable = elem_i.GetGeometry()[node_iter].FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name));
            #pragma omp atomic
            r_nodal_variable[0] += gradient_contribution[0];
            #pragma omp atomic
            r_nodal_variable[1] += gradient_contribution[1];
            #pragma omp atomic
            r_nodal_variable[2] += gradient_contribution[2];                        
            // node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name)) += gradient_contribution;
            node_iter++;
        }

    };

    void CalculateElementThicknessGradients(Element& elem_i, std::string thickness_gradien_name, const ProcessInfo &rCurrentProcessInfo){
        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        Vector u;
        Vector lambda;
        
        // Get state solution
        const auto& rConstElemRef = elem_i;
        rConstElemRef.GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        lambda = 0.5*u;

        Vector d_RHS_d_T;
        double curr_t = elem_i.GetProperties().GetValue(THICKNESS);
        elem_i.GetProperties().SetValue(THICKNESS,1.0);
        elem_i.CalculateRightHandSide(d_RHS_d_T,rCurrentProcessInfo);
        elem_i.GetProperties().SetValue(THICKNESS,curr_t);


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_ppt_d_ft = r_this_geometry[i_node].FastGetSolutionStepValue(D_PPT_D_FT);
            #pragma omp atomic
            r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(thickness_gradien_name)) += d_ppt_d_ft * inner_prod(d_RHS_d_T,lambda) / number_of_nodes;
        }    
    };    

    void CalculateElementObjShapeGradients(Element& elem_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();


        Element::NodesArrayType node_array;
        for (auto& node_i : elem_i.GetGeometry()){
            Element::NodeType::Pointer node_p = node_i.Clone();
            node_array.push_back(node_p);
        }

        Element::Pointer p_elem = elem_i.Create(elem_i.Id(),node_array, elem_i.pGetProperties());
        p_elem->SetData(elem_i.GetData());
        p_elem->Set(Flags(elem_i));
        p_elem->Initialize(rCurrentProcessInfo); 

        double elem_before_pert_stress = CalculateElementValue(*p_elem,rCurrentProcessInfo);

        int node_iter = 0;
        for (auto& node_i : p_elem->GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            double pert_elem_stress = 0;
            
            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            pert_elem_stress = CalculateElementValue(*p_elem,rCurrentProcessInfo);
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = (pert_elem_stress-elem_before_pert_stress) / mDelta;


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            pert_elem_stress = CalculateElementValue(*p_elem,rCurrentProcessInfo);
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = (pert_elem_stress-elem_before_pert_stress) / mDelta;

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            pert_elem_stress = CalculateElementValue(*p_elem,rCurrentProcessInfo);      
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = (pert_elem_stress-elem_before_pert_stress) / mDelta;

            // Assemble sensitivity to node
            array_3d& r_nodal_variable = elem_i.GetGeometry()[node_iter].FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name));
            #pragma omp atomic
            r_nodal_variable[0] += gradient_contribution[0];
            #pragma omp atomic
            r_nodal_variable[1] += gradient_contribution[1];
            #pragma omp atomic
            r_nodal_variable[2] += gradient_contribution[2];                        
            // node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name)) += gradient_contribution;
            node_iter++;
        }

    };

    void CalculateElementObjThicknessGradients(Element& elem_i, std::string thickness_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t number_of_nodes = r_this_geometry.size();

        double curr_thickness = elem_i.GetProperties().GetValue(PT);
        elem_i.GetProperties().SetValue(PT,1.0);
        double elem_thick_grad = CalculateElementValue(elem_i,rCurrentProcessInfo);
        elem_i.GetProperties().SetValue(PT,curr_thickness);

        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_ppt_d_ft = r_this_geometry[i_node].FastGetSolutionStepValue(D_PPT_D_FT);
            #pragma omp atomic
            r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(thickness_gradien_name)) += d_ppt_d_ft * elem_thick_grad / number_of_nodes;
        }

    };    

    void CalculateElementObjMaterialGradients(Element& elem_i, std::string material_gradien_name, const ProcessInfo &rCurrentProcessInfo){
        KRATOS_TRY;
        KRATOS_ERROR << "StressOptResponse: CalculateElementObjMaterialGradients NOT implemented.";
        KRATOS_CATCH("");
    };    

    void CalculateConditionShapeGradients(Condition& cond_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){
        KRATOS_TRY;
        KRATOS_ERROR << "StressOptResponse: CalculateConditionShapeGradients NOT implemented.";
        KRATOS_CATCH("");
    };

    void CalculateElementMaterialGradients(Element& elem_i, std::string material_gradien_name, const ProcessInfo &rCurrentProcessInfo){
        KRATOS_TRY;
        KRATOS_ERROR << "StressOptResponse: CalculateElementMaterialGradients NOT implemented.";
        KRATOS_CATCH("")
    };    

    void GetElementAdjointVector(Element& elem_i, Vector &rAdjoints){

        const GeometryType &rgeom = elem_i.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = elem_i.GetGeometry().WorkingSpaceDimension();
        unsigned int local_size;
        if(mIfShell)
            local_size = num_nodes * 2 * dimension;
        else
            local_size = num_nodes * dimension;

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
            if(mIfShell){
                rAdjoints[index++] =
                    rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_X"));
                rAdjoints[index++] =
                    rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Y"));
                rAdjoints[index++] =
                    rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_ROTATION_Z"));                
            }
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

    LinearSolverType::Pointer mpLinearSystemSolver;
    StrategyType* mpStrategy;
    ModelPart* mpADJModelPart;

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
    bool mIfShell;
    

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
