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

#ifndef PARTITION_MASS_OPT_RESPONSE_H
#define PARTITION_MASS_OPT_RESPONSE_H

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

class PartitionMassOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PartitionMassOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(PartitionMassOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PartitionMassOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings )
        : Response(ResponseName,"interface",rModel, ResponseSettings){
        }

    /// Destructor.
    virtual ~PartitionMassOptResponse()
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
            <<"PartitionMassOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"PartitionMassOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;

        }
    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {

        partitions.clear();
        double total_val = 0.0;
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
			for (auto& elem_i : r_eval_object.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){

                    double element_density = elem_i.GetProperties().GetValue(DENSITY);
                    int element_density_f = 1;
                    int element_density_c = 1;

                    if((element_density>=(element_density_c-0.01))&&(element_density<=element_density_c)){
                        if(partitions.count(element_density_c)<1){
                            partitions.insert({element_density_c,elem_i.GetGeometry().Volume()*element_density});
                        }
                        else{
                            partitions.at(element_density_c) += elem_i.GetGeometry().Volume()*element_density;
                        }
                        total_val += elem_i.GetGeometry().Volume()*element_density;
                    }else if((element_density>=(element_density_f))&&(element_density<=element_density_f+0.01)){
                        if(partitions.count(element_density_f)<1){
                            partitions.insert({element_density_f,elem_i.GetGeometry().Volume()*element_density});
                        }
                        else{
                            partitions.at(element_density_f) += elem_i.GetGeometry().Volume()*element_density;
                        }
                        total_val += elem_i.GetGeometry().Volume()*element_density;
                    }
                }
			}
        }
        // for (auto const &pair: partitions) {
        //     std::cout << " {" << pair.first << ": " << pair.second << "}\n";
        // }
        return total_val;
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
        double element_density_f = std::floor(element_density);
        double element_val =  elem_i.GetGeometry().Volume() * std::pow(element_density-element_density_f,3);

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
                VariableUtils().SetHistoricalVariableToZero(D_PARTITION_MASS_D_FD, controlled_model_part.Nodes());

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

        double elem_dens_grad = 0.0;
        double element_density = elem_i.GetProperties().GetValue(DENSITY);
        int element_density_f = 8;
        int element_density_c = 8;

        if((element_density>=(element_density_c-0.01))&&(element_density<=element_density_c)){
            elem_dens_grad = elem_i.GetGeometry().Volume();
        }else if((element_density>=(element_density_f))&&(element_density<=element_density_f+0.01)){
            elem_dens_grad = elem_i.GetGeometry().Volume();
        }


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pd_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PD_D_FD);
            r_this_geometry[i_node].FastGetSolutionStepValue(D_PARTITION_MASS_D_FD) += d_pd_d_fd * elem_dens_grad / number_of_nodes;
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
        return "PartitionMassOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PartitionMassOptResponse";
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
    std::map<int,double> partitions;

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
//      PartitionMassOptResponse& operator=(PartitionMassOptResponse const& rOther);

    /// Copy constructor.
//      PartitionMassOptResponse(PartitionMassOptResponse const& rOther);


    ///@}

}; // Class PartitionMassOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // PARTITION_MASS_OPT_RESPONSE_H
