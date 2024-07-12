//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

# pragma once

/* System includes */


/* External includes */

/* Project includes */
#include "includes/table_accessor.h"
#include "includes/properties.h"
#include "includes/kratos_parameters.h"


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

///@}
///@name Kratos Classes
///@{

/**
 * @class ReadAndSetAccessorsUtilities
 * @ingroup KratosCore
 * @details static methods class to read and add Accessors to the ReadMaterialsUtility
 * @author Alejandro Cornejo
 */
class KRATOS_API(KRATOS_CORE) ReadAndSetAccessorsUtilities
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void ReadAndSetAccessors(
        const Parameters MaterialData,
        Properties &rProperty);

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
}; /* Class ReadAndSetAccessorsUtilities */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/
