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

#if !defined(KRATOS_ALM_FAST_INIT_PROCESS)
#define KRATOS_ALM_FAST_INIT_PROCESS

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
// This process initializes the variables related with the ALM 
/** Detail class definition.
*/
class ALMFastInit
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMFastInit
    KRATOS_CLASS_POINTER_DEFINITION(ALMFastInit);
    
    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ALMFastInit( ModelPart& rThisModelPart):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY;
        
        KRATOS_CATCH(""); 
    }

    /// Destructor.
    ~ALMFastInit() override
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
        const array_1d<double, 3> zero_vector(3, 0.0);
        
        // We differentiate between frictional or frictionless
        const bool is_frictional = mrThisModelPart.Is(SLIP);
        
        // We initialize the penalty parameter
        const double& epsilon = mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY];
        
        bool init_delta_normal = false;
        Matrix zero_delta_normal;
        if (mrThisModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            init_delta_normal = true;
            const unsigned int dimension = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
            zero_delta_normal = ZeroMatrix( dimension, dimension );
        }
        
        // We iterate over the node
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        
        #pragma omp parallel for firstprivate(zero_vector)
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            // Weighted values
            it_node->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
            if (is_frictional == true)
            {
                it_node->FastGetSolutionStepValue(WEIGHTED_SLIP) = 0.0;
            }
            
            // Penalty parameter
            it_node->SetValue(INITIAL_PENALTY, epsilon);
            
            // Nodal area
            it_node->SetValue(NODAL_AREA, 0.0);
            
            // Auxiliar values
            it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, 0.0);
            if (is_frictional == true)
            {
                it_node->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, 0.0);
            }
            
            // The normal and tangents vectors
            it_node->SetValue(NORMAL, zero_vector);
            
            // The delta normal if necessary
            if (init_delta_normal == true)
            {
                it_node->SetValue(DELTA_NORMAL, zero_delta_normal);
            }
        }
        
        // Now we iterate over the conditions
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for firstprivate(zero_vector)
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            // The normal and tangents vectors
            it_cond->SetValue(NORMAL, zero_vector);
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
        return "ALMFastInit";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALMFastInit";
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
    ALMFastInit& operator=(ALMFastInit const& rOther) = delete;

    /// Copy constructor.
    //ALMFastInit(ALMFastInit const& rOther);


    ///@}

}; // Class ALMFastInit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ALMFastInit& rThis);
// 
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ALMFastInit& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }

}
#endif /* KRATOS_ALM_FAST_INIT_PROCESS defined */
