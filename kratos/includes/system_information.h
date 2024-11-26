//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SystemInformation
 * @ingroup KratosCore
 * @brief This class provides system information
 * @details The implementations come from CppBenchmark, with MIT license, see https://github.com/chronoxor/CppBenchmark/
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) SystemInformation 
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SystemInformation
    KRATOS_CLASS_POINTER_DEFINITION(SystemInformation);

    ///@}
    ///@name Life Cycle
    ///@{
    SystemInformation()
    {}

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SystemInformation" ;

        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "SystemInformation";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {rOStream << "SystemInformation class";}

    ///@}
}; // class
///@}

///@} addtogroup block

} // namespace Kratos