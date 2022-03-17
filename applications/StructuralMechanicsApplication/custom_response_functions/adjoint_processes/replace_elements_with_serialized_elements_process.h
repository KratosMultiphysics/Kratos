// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Mahmoud Sesa
//

#if !defined(KRATOS_REPLACE_ELEMENTS_WITH_SERIALIZED_ELEMENTS_PROCESS)
#define KRATOS_REPLACE_ELEMENTS_WITH_SERIALIZED_ELEMENTS_PROCESS

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ReplaceElementsWithSerializedElementsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReplaceElementsWithSerializedElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceElementsWithSerializedElementsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor.
    ReplaceElementsWithSerializedElementsProcess(
        ModelPart& rMainModelPart, ModelPart& rLoadedModelPart
        ):mrMainModelPart(rMainModelPart), mrLoadedModelPart(rLoadedModelPart)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.

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
        return "ReplaceElementsWithSerializedElementsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ReplaceElementsWithSerializedElementsProcess";
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

    ModelPart& mrMainModelPart;              // The main model part
    ModelPart& mrLoadedModelPart;              // The loaded model part

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
    ReplaceElementsWithSerializedElementsProcess& operator=(ReplaceElementsWithSerializedElementsProcess const& rOther) = delete;

    /// Copy constructor.
    ReplaceElementsWithSerializedElementsProcess(ReplaceElementsWithSerializedElementsProcess const& rOther) = delete;


    ///@}

}; // Class ReplaceElementsWithSerializedElementsProcess

///@}

///@name Type Definitions
///@{


///@}



}
#endif /* KRATOS_REPLACE_ELEMENTS_WITH_SERIALIZED_ELEMENTS_PROCESS defined */
