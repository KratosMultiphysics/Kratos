//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos {

///@addtogroup MetisApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(METIS_APPLICATION) LegacyPartitioningUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LegacyPartitioningUtilities
    KRATOS_CLASS_POINTER_DEFINITION(LegacyPartitioningUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LegacyPartitioningUtilities() = delete;

    /// Copy constructor.
    LegacyPartitioningUtilities(LegacyPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LegacyPartitioningUtilities& operator=(LegacyPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{


    ///@}

}; // Class LegacyPartitioningUtilities

///@}

///@} addtogroup block

} // namespace Kratos
