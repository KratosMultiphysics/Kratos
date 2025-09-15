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

#ifndef LINEAR_STRAIN_ENERGY_OPT_RESPONSE_H
#define LINEAR_STRAIN_ENERGY_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------


// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_responses/response.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"

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

class LinearStrainEnergyOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearStrainEnergyOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(LinearStrainEnergyOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearStrainEnergyOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings )
        : Response(ResponseName,"linear_strain_energy",rModel, ResponseSettings){
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
                        KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;
                }
            }
        }

    /// Destructor.
    virtual ~LinearStrainEnergyOptResponse()
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
            <<"LinearStrainEnergyOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"LinearStrainEnergyOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;

        }
    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        double total_strain_energy = 0.0;
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const ProcessInfo &CurrentProcessInfo = r_eval_object.GetProcessInfo();
            // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
            total_strain_energy += block_for_each<SumReduction<double>>(r_eval_object.Elements(), [&](auto& elem_i) {
                const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
                if(element_is_active)
                {
                    Matrix LHS;
                    Vector RHS;
                    Vector u;

                    // Get state solution relevant for energy calculation
                    const auto& rConstElemRef = elem_i;
                    rConstElemRef.GetValuesVector(u,0);

                    elem_i.CalculateLocalSystem(LHS,RHS,CurrentProcessInfo);

                    // Compute strain energy
                    return 0.5 * inner_prod(u,prod(LHS,u));
                } else {
                    return 0.0;
                }
            });
        }
        return total_strain_energy;
    };


    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        for(long unsigned int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            const ProcessInfo &CurrentProcessInfo = controlled_model_part.GetProcessInfo();
            auto control_type = mrResponseSettings["control_types"][i].GetString();
            std::string grad_field_name;
            if(control_type=="shape"){
                grad_field_name = mrResponseSettings["gradient_settings"]["shape_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<array_1d<double,3>>>::Get(grad_field_name), controlled_model_part.Nodes());
                block_for_each(controlled_model_part.Elements(), [&](auto& elem_i) {
                    if(elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true)
                        CalculateElementShapeGradients(elem_i,grad_field_name,CurrentProcessInfo);
                });

                for (auto& cond_i : controlled_model_part.Conditions())
                    if(cond_i.IsDefined(ACTIVE) ? cond_i.Is(ACTIVE) : true)
                        CalculateConditionShapeGradients(cond_i,grad_field_name,CurrentProcessInfo);
            }
            else if(control_type=="material"){
                grad_field_name = mrResponseSettings["gradient_settings"]["material_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<double>>::Get(grad_field_name), controlled_model_part.Nodes());
                block_for_each(controlled_model_part.Elements(), [&](auto& elem_i) {
                    if(elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true)
                        CalculateElementMaterialGradients(elem_i,grad_field_name,CurrentProcessInfo);
                });
            }
            else if(control_type=="thickness"){
                grad_field_name = mrResponseSettings["gradient_settings"]["thickness_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<double>>::Get(grad_field_name), controlled_model_part.Nodes());
                block_for_each(controlled_model_part.Elements(), [&](auto& elem_i) {
                    if(elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true)
                        CalculateElementThicknessGradients(elem_i,grad_field_name,CurrentProcessInfo);
                });
            }

        }

		KRATOS_CATCH("");

    };

    void CalculateElementShapeGradients(Element& elem_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        Element::NodesArrayType node_array;
        for (auto& node_i : elem_i.GetGeometry()){
            Element::NodeType::Pointer node_p = node_i.Clone();
            node_array.push_back(node_p);
        }

        Element::Pointer p_elem = elem_i.Create(elem_i.Id(),node_array, elem_i.pGetProperties());
        p_elem->SetData(elem_i.GetData());
        p_elem->Set(Flags(elem_i));
        p_elem->Initialize(rCurrentProcessInfo);

        Vector u;
        Vector lambda;
        Vector RHS;

        // Get state solution
        p_elem->GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        lambda = 0.5*u;

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
            AtomicAdd(r_nodal_variable[0], gradient_contribution[0]);
            AtomicAdd(r_nodal_variable[1], gradient_contribution[1]);
            AtomicAdd(r_nodal_variable[2], gradient_contribution[2]);
            // node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name)) += gradient_contribution;
            node_iter++;
        }

    };

    void CalculateConditionShapeGradients(Condition& cond_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        Element::NodesArrayType node_array;
        for (auto& node_i : cond_i.GetGeometry()){
            Element::NodeType::Pointer node_p = node_i.Clone();
            node_array.push_back(node_p);
        }

        Condition::Pointer p_cond = cond_i.Create(cond_i.Id(),node_array, cond_i.pGetProperties());
        p_cond->SetData(cond_i.GetData());
        p_cond->Set(Flags(cond_i));
        p_cond->Initialize(rCurrentProcessInfo);

        Vector u;
        Vector lambda;
        Vector RHS;

        // Get state solution
        p_cond->GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        lambda = 0.5*u;

        // Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
        p_cond->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        int node_iter = 0;
        for (auto& node_i : p_cond->GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            Vector RHS_perturbed = Vector(RHS.size());
            Vector derived_RHS = Vector(RHS.size());

            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            p_cond->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = inner_prod(lambda, derived_RHS);


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            p_cond->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = inner_prod(lambda, derived_RHS);

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            p_cond->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = inner_prod(lambda, derived_RHS);

            // Assemble sensitivity to node
            array_3d& r_nodal_variable = cond_i.GetGeometry()[node_iter].FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name));
            AtomicAdd(r_nodal_variable[0], gradient_contribution[0]);
            AtomicAdd(r_nodal_variable[1], gradient_contribution[1]);
            AtomicAdd(r_nodal_variable[2], gradient_contribution[2]);
            // node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name)) += gradient_contribution;
            node_iter++;
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
        lambda = 0.5*u;

        Vector d_RHS_d_E;
        double curr_e = elem_i.GetProperties().GetValue(YOUNG_MODULUS);
        elem_i.GetProperties().SetValue(YOUNG_MODULUS,1.0);
        elem_i.CalculateRightHandSide(d_RHS_d_E,rCurrentProcessInfo);
        elem_i.GetProperties().SetValue(YOUNG_MODULUS,curr_e);


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pe_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PE_D_FD);
            AtomicAdd(r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(material_gradien_name)), d_pe_d_fd * inner_prod(d_RHS_d_E,lambda) / number_of_nodes);
        }

    };


    void CalculateElementThicknessGradients(Element& elem_i, std::string thickness_gradien_name, const ProcessInfo &rCurrentProcessInfo){
        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
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
            AtomicAdd(r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(thickness_gradien_name)), d_ppt_d_ft * inner_prod(d_RHS_d_T,lambda) / number_of_nodes);
        }
    };

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
        return "LinearStrainEnergyOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LinearStrainEnergyOptResponse";
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
//      LinearStrainEnergyOptResponse& operator=(LinearStrainEnergyOptResponse const& rOther);

    /// Copy constructor.
//      LinearStrainEnergyOptResponse(LinearStrainEnergyOptResponse const& rOther);


    ///@}

}; // Class LinearStrainEnergyOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // LINEAR_STRAIN_ENERGY_OPT_RESPONSE_H
