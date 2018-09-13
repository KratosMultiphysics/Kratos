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

/**
 * @class ALMFastInit
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process initializes the variables related with the ALM
 * @details Initializes:
 * - SLIP flag
 * - Nodal INITIAL_PENALTY
 * - Nodal WEIGHTED_GAP
 * - Nodal WEIGHTED_SLIP
 * - Nodal DYNAMIC_FACTOR
 * - Nodal AUGMENTED_NORMAL_CONTACT_PRESSURE
 * - Nodal AUGMENTED_TANGENT_CONTACT_PRESSURE
 * - Condition NORMAL
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ALMFastInit
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

    void Execute() override;

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

    ModelPart& mrThisModelPart; /// The model part to initialize

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
