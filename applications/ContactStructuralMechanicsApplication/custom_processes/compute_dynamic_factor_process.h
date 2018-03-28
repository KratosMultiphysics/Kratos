
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

#if !defined(KRATOS_COMPUTE_DYNAMIC_FACTOR_PROCESS)
#define KRATOS_COMPUTE_DYNAMIC_FACTOR_PROCESS

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
 * @class ComputeDynamicFactorProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process is used in order to compute the dynamic factor for dynamic problems
 * @details The estimation is done in proportion to the weighthed gap (current and previous)
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputeDynamicFactorProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeDynamicFactorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeDynamicFactorProcess);

    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeDynamicFactorProcess( ModelPart& rThisModelPart):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~ComputeDynamicFactorProcess() override
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
        return "ComputeDynamicFactorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeDynamicFactorProcess";
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
    ComputeDynamicFactorProcess& operator=(ComputeDynamicFactorProcess const& rOther) = delete;

    /// Copy constructor.
    //ComputeDynamicFactorProcess(ComputeDynamicFactorProcess const& rOther);


    ///@}

}; // Class ComputeDynamicFactorProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ComputeDynamicFactorProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ComputeDynamicFactorProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_COMPUTE_DYNAMIC_FACTOR_PROCESS defined */
