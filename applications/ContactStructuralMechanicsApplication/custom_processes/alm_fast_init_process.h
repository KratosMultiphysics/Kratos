// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
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
    virtual ~ALMFastInit()
    {
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
    
    virtual void Execute() override
    {
        KRATOS_TRY;
        
        // We initialize the zero vector
        const array_1d<double, 3> zerovector(3, 0.0);
        
        // We differentiate between frictional or frictionless
        const bool IsFrictional = mrThisModelPart.Is(SLIP);
        
        bool InitDeltaNormal = false;
        Matrix ZeroDeltaNormal;
        if (mrThisModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            InitDeltaNormal = true;
            const unsigned int Dimension = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
            ZeroDeltaNormal = ZeroMatrix( Dimension, Dimension );
        }
        
        // We iterate over the node
        NodesArrayType& pNodes = mrThisModelPart.Nodes();
        auto numNodes = pNodes.end() - pNodes.begin();
        
        #pragma omp parallel for firstprivate(zerovector)
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNodes.begin() + i;
            
            // Weighted values
            itNode->FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
            if (IsFrictional == true)
            {
                itNode->FastGetSolutionStepValue(WEIGHTED_SLIP) = 0.0;
            }
            
            // Nodal area
            itNode->SetValue(NODAL_AREA, 0.0);
            
            // Auxiliar values
            itNode->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, 0.0);
            if (IsFrictional == true)
            {
                itNode->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, 0.0);
            }
            
            // The normal and tangents vectors
            itNode->SetValue(NORMAL, zerovector);
            
            // The delta normal if necessary
            if (InitDeltaNormal == true)
            {
                itNode->SetValue(DELTA_NORMAL, ZeroDeltaNormal);
            }
            
            
        }
        
        // Now we iterate over the conditions
        ConditionsArrayType& pConditions = mrThisModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();
        
        #pragma omp parallel for firstprivate(zerovector)
        for(int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            // The normal and tangents vectors
            itCond->SetValue(NORMAL, zerovector);
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
    virtual std::string Info() const override
    {
        return "ALMFastInit";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALMFastInit";
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
    ALMFastInit& operator=(ALMFastInit const& rOther);

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
