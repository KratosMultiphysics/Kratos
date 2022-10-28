// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef PARTITION_INTERFACE_STRESS_OPT_RESPONSE_H
#define PARTITION_INTERFACE_STRESS_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_responses/response.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/variable_utils.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "processes/find_conditions_neighbours_process.h"
#include "response_functions/adjoint_response_function.h"
#include "custom_elements/adjoint_small_displacement_element.h"
#include "custom_strategies/strategies/helmholtz_strategy.h"

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

class KRATOS_API(OPTIMIZATION_APPLICATION) PartitionInterfaceStressOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef Variable<double> array_1d_component_type;  
    typedef array_1d<double,3> array_3d;
    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;    
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;    
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;  
    typedef HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> StrategyType;    

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
        for(int i=0;i<mrResponseSettings["evaluated_objects"].size();i++){
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
            const std::size_t domain_size = r_eval_object.GetProcessInfo()[DOMAIN_SIZE];
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

    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        // compute sensitivities       
        for(int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            const ProcessInfo &CurrentProcessInfo = controlled_model_part.GetProcessInfo();
            auto control_type = mrResponseSettings["control_types"][i].GetString();
            VariableUtils().SetHistoricalVariableToZero(D_STRESS_D_FD, controlled_model_part.Nodes());

            //obj
			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    CalculateElementObjMaterialGradients(elem_i,CurrentProcessInfo);                                               
                }
            }
        }

		KRATOS_CATCH("");
 
    };

    void CalculateElementObjMaterialGradients(Element& elem_i, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();
      
        std::vector<double> gp_weights_vector;
        elem_i.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, gp_weights_vector, rCurrentProcessInfo);        
        std::vector<double> stress_gp_vector;
        elem_i.CalculateOnIntegrationPoints(KratosComponents<Variable<double>>::Get("VON_MISES_STRESS"), stress_gp_vector, rCurrentProcessInfo); 
        double element_density = elem_i.GetProperties().GetValue(DENSITY);

        double elem_sens = 0.0;
        for(IndexType i = 0; i < gp_weights_vector.size(); i++){            
            elem_sens += gp_weights_vector[i] * stress_gp_vector[i] * 2 * std::pow((element_density) * (1.0-element_density),1) * (1.0 - 2.0 * element_density);
        }

        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pd_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PD_D_FD);
            r_this_geometry[i_node].FastGetSolutionStepValue(D_STRESS_D_FD) += d_pd_d_fd * elem_sens / number_of_nodes;
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
