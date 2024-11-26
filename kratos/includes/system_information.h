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

    /**
     * @brief Get OS version string
     * @return The OS version considered
     */
    static std::string OSVersion();

    /**
     * @brief CPU architecture string
     * @return The CPU architecture of the current machine
     */
    static std::string CPUArchitecture();

    /**
     * @brief CPU logical cores count
     * @return The number of logical CPU cores
     */
    static std::size_t CPULogicalCores();

    /**
     * @brief CPU physical cores count
     * @return The number of physical CPU cores
     */
    static std::size_t CPUPhysicalCores();

    /**
     * @brief CPU clock speed in Hz
     * @return The CPU clock speed
     */
    static std::size_t CPUClockSpeed();

    /**
     * @brief Is CPU Hyper-Threading enabled?
     * @return If the hyperthreading is enabled
     */
    static bool CPUHyperThreading();

    /**
     * @brief Total RAM in bytes
     * @return The total RAM in bytes
     */
    static std::size_t RamTotal();

    /**
     * @brief Free RAM in bytes
     * @return The free RAM in bytes
     */
    static std::size_t RamFree();

    /**
     * @brief Generate clock speed string
     * @details Will return a pretty string of Hz, kHz, MHz, GHz based on the given clock speed in hertz.
     * @param Hertz Clock speed value in hertz
     * @return String with clock speed representation
     */
    static std::string GenerateClockSpeed(const std::size_t Hertz);
    
    /**
     * @brief Generate data size string
     * @details Will return a pretty string of bytes, KiB, MiB, GiB, TiB based on the given bytes.
     * @param Bytes Data size in bytes
     * @return String with data size representation
     */
    static std::string GenerateDataSize(const std::size_t Bytes);

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