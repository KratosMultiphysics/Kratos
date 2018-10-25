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

#if !defined(KRATOS_TOTAL_STRUCTURAL_MASS_PROCESS)
#define KRATOS_TOTAL_STRUCTURAL_MASS_PROCESS

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
 * @class TotalStructuralMassProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes the total mass of a structure
 * @details It takes into account the noda-mass, beam, shells and solid elements
 *
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TotalStructuralMassProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TotalStructuralMassProcess
    KRATOS_CLASS_POINTER_DEFINITION(TotalStructuralMassProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TotalStructuralMassProcess(
        ModelPart& rThisModelPart
        ):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~TotalStructuralMassProcess() override
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

    static double CalculateElementMass(Element& rElement, const std::size_t DomainSize);

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
        return "TotalStructuralMassProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TotalStructuralMassProcess";
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
    TotalStructuralMassProcess& operator=(TotalStructuralMassProcess const& rOther) = delete;

    /// Copy constructor.
    //TotalStructuralMassProcess(TotalStructuralMassProcess const& rOther);


    ///@}

}; // Class TotalStructuralMassProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   TotalStructuralMassProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const TotalStructuralMassProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_TOTAL_STRUCTURAL_MASS_PROCESS defined */
