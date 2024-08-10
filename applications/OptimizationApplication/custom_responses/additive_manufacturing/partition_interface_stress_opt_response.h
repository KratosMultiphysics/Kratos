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

#ifndef PARTITION_INTERFACE_STRESS_OPT_RESPONSE_H
#define PARTITION_INTERFACE_STRESS_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------


// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_responses/response.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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

class PartitionInterfaceStressOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PartitionInterfaceStressOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(PartitionInterfaceStressOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PartitionInterfaceStressOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings)
        : Response(ResponseName,"stress",rModel, ResponseSettings){

        }

    /// Destructor.
    virtual ~PartitionInterfaceStressOptResponse()
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
            <<"PartitionInterfaceStressOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"PartitionInterfaceStressOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;

        }
    };

    // --------------------------------------------------------------------------
    double CalculateElementStress(Element& elem_i, const ProcessInfo &rCurrentProcessInfo){

        std::vector<double> gp_weights_vector;
        std::vector<double> stress_gp_vector;
        elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, rCurrentProcessInfo);
        elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), stress_gp_vector, rCurrentProcessInfo);
        double element_density = elem_i.GetProperties().GetValue(DENSITY);

        double elem_value = 0.0;
        for(IndexType i = 0; i < gp_weights_vector.size(); i++){
            elem_value += gp_weights_vector[i] * stress_gp_vector[i] * std::pow((element_density) * (1.0-element_density),2);
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
            intg_stress += block_for_each<SumReduction<double>>(r_eval_object.Elements(), [&](auto& elem_i) {
                const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
                if(element_is_active){
                    return CalculateElementStress(elem_i,CurrentProcessInfo);
                } else {
                    return 0.0;
                }
            });
        }
        return intg_stress;
    };

    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        // compute sensitivities
        for(long unsigned int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            const ProcessInfo &CurrentProcessInfo = controlled_model_part.GetProcessInfo();
            VariableUtils().SetHistoricalVariableToZero(D_STRESS_D_FD, controlled_model_part.Nodes());

            block_for_each(controlled_model_part.Elements(), [&](auto& elem_i) {
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    auto& r_this_geometry = elem_i.GetGeometry();
                    const std::size_t number_of_nodes = r_this_geometry.size();

                    std::vector<double> gp_weights_vector;
                    elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, CurrentProcessInfo);
                    std::vector<double> stress_gp_vector;
                    elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), stress_gp_vector, CurrentProcessInfo);
                    double element_density = elem_i.GetProperties().GetValue(DENSITY);

                    double elem_sens = 0.0;
                    for(IndexType i = 0; i < gp_weights_vector.size(); i++){
                        elem_sens += gp_weights_vector[i] * stress_gp_vector[i] * 2 * std::pow((element_density) * (1.0-element_density),1) * (1.0 - 2.0 * element_density);
                    }

                    for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
                        const auto& d_pd_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PD_D_FD);
                        AtomicAdd(r_this_geometry[i_node].FastGetSolutionStepValue(D_STRESS_D_FD), d_pd_d_fd * elem_sens / number_of_nodes);
                    }
                }
            });
        }

		KRATOS_CATCH("");

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
        return "PartitionInterfaceStressOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PartitionInterfaceStressOptResponse";
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
//      PartitionInterfaceStressOptResponse& operator=(PartitionInterfaceStressOptResponse const& rOther);

    /// Copy constructor.
//      PartitionInterfaceStressOptResponse(PartitionInterfaceStressOptResponse const& rOther);


    ///@}

}; // Class PartitionInterfaceStressOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // PARTITION_INTERFACE_STRESS_OPT_RESPONSE_H
