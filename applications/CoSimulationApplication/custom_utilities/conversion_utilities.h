//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

#if !defined(KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED )
#define  KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CO_SIMULATION_APPLICATION) ConversionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConversionUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ConversionUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConversionUtilities() = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void ConvertElementalDataToNodalData(ModelPart& rModelPart);

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ConversionUtilities& operator=(ConversionUtilities const& rOther);

    /// Copy constructor.
    ConversionUtilities(ConversionUtilities const& rOther);


    ///@}

}; // Class ConversionUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED  defined
