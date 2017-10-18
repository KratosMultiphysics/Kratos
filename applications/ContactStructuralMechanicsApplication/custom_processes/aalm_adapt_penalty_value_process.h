 
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

#if !defined(KRATOS_AALM_ADAPT_PENALTY_VALUE_PROCESS)
#define KRATOS_AALM_ADAPT_PENALTY_VALUE_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"

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
// This process adapts the penalty following the algorithm (Algorithm 3) from "The adapted augmented Lagrangian method: a new method for the resolution of the mechanical frictional contact problem" Philippe Bussetta · Daniel Marceau ·Jean-Philippe Ponthot
/** Detail class definition.
*/
class AALMAdaptPenaltyValueProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AALMAdaptPenaltyValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(AALMAdaptPenaltyValueProcess);
    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AALMAdaptPenaltyValueProcess( ModelPart& rThisModelPart):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY;
        
        KRATOS_CATCH(""); 
    }

    /// Destructor.
    ~AALMAdaptPenaltyValueProcess() override
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
        KRATOS_TRY;
        
        // We initialize the zero vector
        const double& penalty_parameter = mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY];
        const double& max_gap_factor = mrThisModelPart.GetProcessInfo()[MAX_GAP_FACTOR];
        
        // We iterate over the node
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            // Weighted values
            const double& current_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP);
            const double& previous_gap = it_node->FastGetSolutionStepValue(WEIGHTED_GAP, 1);
            
            // Nodal H
            const double& nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);
            const double max_gap = max_gap_factor * nodal_h; // NOTE: This value must be studied
            
            if ((current_gap * previous_gap) < 0.0)
            {
                if (previous_gap > max_gap)
                {
                    it_node->SetValue(INITIAL_PENALTY, std::abs(penalty_parameter * previous_gap / (current_gap) * (std::abs(current_gap) + max_gap)/(current_gap - previous_gap)));
                }
                else
                {
                    it_node->SetValue(INITIAL_PENALTY, std::abs(penalty_parameter * previous_gap / (10.0 * current_gap)));
                }
            }
            else if (current_gap > max_gap)
            {
                if (std::abs(current_gap - previous_gap) > std::max(current_gap/10.0, std::max(previous_gap/1.0, 5 * max_gap)))
                {
                    it_node->SetValue(INITIAL_PENALTY, 2.0 * penalty_parameter);
                }
                else if ((std::abs(current_gap) <= std::abs(previous_gap) * 1.01 || std::abs(current_gap) >= std::abs(previous_gap) * 0.99) && (std::abs(current_gap) < 10.0 *  max_gap))
                {
                    it_node->SetValue(INITIAL_PENALTY, penalty_parameter * std::pow((std::sqrt(std::abs(current_gap)/max_gap - 1.0) + 1.0), 2.0));
                }
                else if (std::abs(current_gap) > std::abs(previous_gap) * 1.01)
                {
                    it_node->SetValue(INITIAL_PENALTY, 2.0 * penalty_parameter * (previous_gap/current_gap));
                }
                else
                {
                    it_node->SetValue(INITIAL_PENALTY, penalty_parameter * (std::sqrt(std::abs(current_gap)/max_gap - 1.0) + 1.0));
                }
            }
            else
            {
                it_node->SetValue(INITIAL_PENALTY, penalty_parameter);
            }
        }

        KRATOS_CATCH("");
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
        return "AALMAdaptPenaltyValueProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AALMAdaptPenaltyValueProcess";
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
    
    ModelPart& mrThisModelPart;

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
    AALMAdaptPenaltyValueProcess& operator=(AALMAdaptPenaltyValueProcess const& rOther) = delete;

    /// Copy constructor.
    //AALMAdaptPenaltyValueProcess(AALMAdaptPenaltyValueProcess const& rOther);


    ///@}

}; // Class AALMAdaptPenaltyValueProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   AALMAdaptPenaltyValueProcess& rThis);
// 
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const AALMAdaptPenaltyValueProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }

}
#endif /* KRATOS_AALM_ADAPT_PENALTY_VALUE_PROCESS defined */
