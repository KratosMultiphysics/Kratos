// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes


// External includes


// Project includes
#include "includes/model_part.h"

// Application includes


namespace Kratos::Testing {

///@addtogroup ConvectionDiffusionApplication
///@{

///@name Kratos Classes
///@{

/// Conveciton diffusion testing utilities.
/** This class contains the auxiliary functions used in the application tests.
 * Note that this class must never be included outside the Kratos::Testing namespace.
*/
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) ConvectionDiffusionTestingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConvectionDiffusionTestingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ConvectionDiffusionTestingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConvectionDiffusionTestingUtilities() = delete;

    /// Copy constructor.
    ConvectionDiffusionTestingUtilities(ConvectionDiffusionTestingUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ConvectionDiffusionTestingUtilities& operator=(ConvectionDiffusionTestingUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static void SetEntityUnitTestModelPart(ModelPart &rModelPart);

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


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

}; // Class ConvectionDiffusionTestingUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                ConvectionDiffusionTestingUtilities& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ConvectionDiffusionTestingUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos::Testing.
