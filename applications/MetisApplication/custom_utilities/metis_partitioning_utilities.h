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

/// This class contains utilities used by the metis partitioners.
/** Detail class definition.
*/
class KRATOS_API(METIS_APPLICATION) MetisPartitioningUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MetisPartitioningUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MetisPartitioningUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MetisPartitioningUtilities() = delete;

    /// Copy constructor.
    MetisPartitioningUtilities(MetisPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MetisPartitioningUtilities& operator=(MetisPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{


    ///@}

}; // Class MetisPartitioningUtilities

///@}

///@} addtogroup block

} // namespace Kratos
