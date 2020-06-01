//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_PROCESS_FACTORY_H_INCLUDED )
#define  KRATOS_PROCESS_FACTORY_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/kratos_components.h"
#include "processes/process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Creates Processes
namespace ProcessFactory
{
    /// Checks if the process is registered
    bool Has(const std::string& ProcessName)
    {
        return KratosComponents< Process >::Has(ProcessName);
    }

    /// Checks if the process is registered
    typename Process::Pointer Create(
        const std::string& ProcessName, Model& rModel, const Parameters ModelParameters)
    {
        KRATOS_ERROR_IF_NOT(Has(ProcessName))
            << "Trying to construct a process: "
            << ProcessName << "\" which does not exist.\n"
            << "The available options (for currently loaded applications) are:\n"
            << KratosComponents< Process >() << std::endl;

        Process const& r_clone_process = KratosComponents< Process >::Get(ProcessName);
        return r_clone_process.Create(rModel, ModelParameters);
    }

}

}  // namespace Kratos.

#endif // KRATOS_PROCESS_FACTORY_H_INCLUDED  defined

