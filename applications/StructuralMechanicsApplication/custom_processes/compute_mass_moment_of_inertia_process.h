// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher
//                   Salman Yousaf
//

#if !defined(KRATOS_COMPUTE_MOMENT_OF_INERTIA_PROCESS)
#define KRATOS_COMPUTE_MOMENT_OF_INERTIA_PROCESS

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
 * @class ComputeMassMomentOfInertiaProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes the moment of inertia (rotational)
 * @details It takes into account all elements in the ModelPart
 *
 * @author Philipp Bucher, Salman Yousaf
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ComputeMassMomentOfInertiaProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeMassMomentOfInertiaProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeMassMomentOfInertiaProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeMassMomentOfInertiaProcess(
        ModelPart& rThisModelPart,
        const Point& rPoint1,
        const Point& rPoint2
        ):mrThisModelPart(rThisModelPart) , mrPoint1(rPoint1), mrPoint2(rPoint2)
    { }

    /// Destructor.
    ~ComputeMassMomentOfInertiaProcess() override = default;

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
        return "ComputeMassMomentOfInertiaProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeMassMomentOfInertiaProcess";
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

    ModelPart& mrThisModelPart; // The main model part
    const Point& mrPoint1;      // Point 1 of the axis of rotation
    const Point& mrPoint2;      // Point 2 of the axis of rotation

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
    ComputeMassMomentOfInertiaProcess& operator=(ComputeMassMomentOfInertiaProcess const& rOther) = delete;

    /// Copy constructor.
    ComputeMassMomentOfInertiaProcess(ComputeMassMomentOfInertiaProcess const& rOther) = delete;


    ///@}

}; // Class ComputeMassMomentOfInertiaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ComputeMassMomentOfInertiaProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ComputeMassMomentOfInertiaProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_COMPUTE_MOMENT_OF_INERTIA_PROCESS defined */
