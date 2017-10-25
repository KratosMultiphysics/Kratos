// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ALM_VARIABLES_CALCULATION_PROCESS)
#define KRATOS_ALM_VARIABLES_CALCULATION_PROCESS

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"
#include "custom_utilities/contact_utilities.h"

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
    
/// Short class definition.
// This process calculates the variables related with the ALM 
/** Detail class definition.
*/
class ALMVariablesCalculationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMVariablesCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(ALMVariablesCalculationProcess);
    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Point                                        PointType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef Geometry<PointType>                     GeometryPointType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ALMVariablesCalculationProcess(
        ModelPart& rThisModelPart,
        Variable<double>& rNodalLengthVariable = NODAL_H,
        Parameters ThisParameters =  Parameters(R"({})")
        ):mrThisModelPart(rThisModelPart), 
          mrNodalLengthVariable(rNodalLengthVariable)
    {
        KRATOS_TRY
        
        Parameters default_parameters = Parameters(R"(
        {
            "stiffness_factor"                     : 10.0,
            "penalty_scale_factor"                 : 1.0
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);
        
        mFactorStiffness = ThisParameters["stiffness_factor"].GetDouble();
        mPenaltyScale = ThisParameters["penalty_scale_factor"].GetDouble();
        
        if (rThisModelPart.GetNodalSolutionStepVariablesList().Has( rNodalLengthVariable ) == false ) // TODO: Ask Riccardo if it is necessary to use GetNodalSolutionStepVariablesList (REQUIRES BUFFER!!!!)
        {
            KRATOS_ERROR << "Missing variable " << rNodalLengthVariable;
        }
        
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ALMVariablesCalculationProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    void Execute() override
    {
        KRATOS_TRY
        
        /* We compute the penalty factor */
        
        // We initialize the mean values
        double mean_young_modulus_slave  = 0.0;
        double mean_nodal_h_slave        = 0.0;
        double mean_young_modulus_master = 0.0;
        double mean_nodal_h_master       = 0.0;
        
        // We initialize the total areas and volumes
        double total_volume_slave  = 0.0;
        double total_area_slave    = 0.0;
        double total_volume_master = 0.0;
        double total_area_master   = 0.0;
        
        // Now we iterate over the conditions to calculate the nodal area
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            // We get the condition geometry
            GeometryType& r_this_geometry = it_cond->GetGeometry();
            const unsigned num_nodes_geometry = r_this_geometry.size();
            
            // We get the values from the element
            Element::Pointer p_elem = it_cond->GetValue(ELEMENT_POINTER);
            
            const double young_modulus = p_elem->GetProperties()[YOUNG_MODULUS];
            const double element_volume = p_elem->GetGeometry().Area();
            
            // We get the values from the condition
            const double condition_area = r_this_geometry.Area();
            const double nodal_condition_area = condition_area/num_nodes_geometry;
            
            if (it_cond->Is(SLAVE) == true)
            {
                #pragma omp atomic
                total_volume_slave += element_volume;
                #pragma omp atomic
                total_area_slave += condition_area;
                #pragma omp atomic
                mean_young_modulus_slave += young_modulus * element_volume;
                
                for (unsigned int i_node = 0; i_node < num_nodes_geometry; i_node++)
                {
                    #pragma omp atomic
                    mean_nodal_h_slave += r_this_geometry[i_node].FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
                }
            }
            
            if (it_cond->Is(MASTER) == true)
            {
                #pragma omp atomic
                total_volume_master += element_volume;
                #pragma omp atomic
                total_area_master += condition_area;
                #pragma omp atomic
                mean_young_modulus_master += young_modulus * element_volume;
                
                for (unsigned int i_node = 0; i_node < num_nodes_geometry; i_node++)
                {
                    #pragma omp atomic
                    mean_nodal_h_master += r_this_geometry[i_node].FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
                }
            }
        }
        
        // Now we divide between the total areas and volumes
        mean_nodal_h_slave /= (total_area_slave + 1.0e-12); 
        mean_young_modulus_slave /= (total_volume_slave + 1.0e-12);
        
        mean_nodal_h_master /= (total_area_master  + 1.0e-12); 
        mean_young_modulus_master /= (total_volume_master + 1.0e-12);
        
        // Finally we compute the penalty factor
        const double penalty_parameter_slave  = mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + 1.0e-12);
        const double scale_factor_slave    = mPenaltyScale * mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + 1.0e-12);
        const double penalty_parameter_master = mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + 1.0e-12);
        const double scale_factor_master   = mPenaltyScale * mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + 1.0e-12); 
        
        mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY] = (penalty_parameter_slave > penalty_parameter_master) ? penalty_parameter_slave : penalty_parameter_master; // NOTE: > or <? , we are supposed to take the largest of the values (more stiff)
        mrThisModelPart.GetProcessInfo()[SCALE_FACTOR] = (scale_factor_slave > scale_factor_master) ? scale_factor_slave : scale_factor_master;
        
        KRATOS_CATCH("")
    }
    
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
    std::string Info() const override
    {
        return "ALMVariablesCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALMVariablesCalculationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    
    ModelPart& mrThisModelPart;              // The main model part
    Variable<double>& mrNodalLengthVariable; // The variable used to messure the lenght of the element
    double mFactorStiffness;                 // The proportion between stiffness and penalty/scale factor
    double mPenaltyScale;                    // The penalty/scale factor proportion

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


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
    ALMVariablesCalculationProcess& operator=(ALMVariablesCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    //ALMVariablesCalculationProcess(ALMVariablesCalculationProcess const& rOther);


    ///@}

}; // Class ALMVariablesCalculationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ALMVariablesCalculationProcess& rThis);
// 
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ALMVariablesCalculationProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }

}
#endif /* KRATOS_ALM_VARIABLES_CALCULATION_PROCESS defined */
