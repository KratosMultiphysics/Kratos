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

#ifndef INTERFACE_OPT_RESPONSE_H
#define INTERFACE_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------


// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_responses/response.h"

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

class InterfaceOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InterfaceOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings )
        : Response(ResponseName,"interface",rModel, ResponseSettings){
        }

    /// Destructor.
    virtual ~InterfaceOptResponse()
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
            <<"InterfaceOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"InterfaceOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;

        }
    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        double total_mass = 0.0;
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const std::size_t domain_size = r_eval_object.GetProcessInfo()[DOMAIN_SIZE];
			for (auto& elem_i : r_eval_object.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active)
                    total_mass += CalculateElementValue(elem_i,domain_size);
			}
        }
        return total_mass;
    };

    double CalculateElementValue(Element& elem_i, const std::size_t DomainSize){
        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t number_of_nodes = r_this_geometry.size();

        // We copy the current coordinates and move the coordinates to the initial configuration
        std::vector<array_1d<double, 3>> current_coordinates(number_of_nodes);
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(current_coordinates[i_node]) = r_this_geometry[i_node].Coordinates();
            noalias(r_this_geometry[i_node].Coordinates()) = r_this_geometry[i_node].GetInitialPosition().Coordinates();
        }

        double element_density = elem_i.GetProperties().GetValue(DENSITY);
        double element_val =  elem_i.GetGeometry().Volume() * std::pow((element_density) * (1.0-element_density),2);

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];
        }

        return element_val;
    }

    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        for(long unsigned int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            const std::size_t domain_size = controlled_model_part.GetProcessInfo()[DOMAIN_SIZE];
            auto control_type = mrResponseSettings["control_types"][i].GetString();
            if(control_type=="material")
                VariableUtils().SetHistoricalVariableToZero(D_INTERFACE_D_FD, controlled_model_part.Nodes());

			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="material")
                        CalculateElementMaterialGradients(elem_i,domain_size);
                }
            }


        }

		KRATOS_CATCH("");

    };

    void CalculateElementMaterialGradients(Element& elem_i, const std::size_t DomainSize){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t number_of_nodes = r_this_geometry.size();

        double curr_density = elem_i.GetProperties().GetValue(DENSITY);
        double elem_dens_grad = elem_i.GetGeometry().Volume() * 2 * std::pow((curr_density) * (1.0-curr_density),1) * (1.0 - 2.0 * curr_density);


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pd_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PD_D_FD);
            r_this_geometry[i_node].FastGetSolutionStepValue(D_INTERFACE_D_FD) += d_pd_d_fd * elem_dens_grad / number_of_nodes;
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
        return "InterfaceOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceOptResponse";
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
//      InterfaceOptResponse& operator=(InterfaceOptResponse const& rOther);

    /// Copy constructor.
//      InterfaceOptResponse(InterfaceOptResponse const& rOther);


    ///@}

}; // Class InterfaceOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // INTERFACE_OPT_RESPONSE_H
