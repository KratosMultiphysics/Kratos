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

// System includes

// External includes

// Project includes
#include "custom_utilities/process_factory_utility.h"

namespace Kratos
{
ProcessFactoryUtility::ProcessFactoryUtility(ListType& ProcessesList)
{
    const std::size_t size_processes = len(ProcessesList);
    mProcesses.resize(size_processes);
    for (std::size_t i_process = 0; i_process < size_processes; ++i_process)
        mProcesses[i_process] = pybind11::cast<ObjectType>(ProcessesList[i_process]);
}

/***********************************************************************************/
/***********************************************************************************/

ProcessFactoryUtility::ProcessFactoryUtility(ObjectType& rProcess)
{
    mProcesses.resize(1);
    mProcesses[0] = rProcess;
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::AddProcess(ObjectType& rProcess)
{
    mProcesses.push_back(rProcess);
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::AddProcesses(ListType& ProcessesList)
{
    const std::size_t size_processes = len(ProcessesList);
    for (std::size_t i_process = 0; i_process < size_processes; ++i_process)
        mProcesses.push_back(pybind11::cast<ObjectType>(ProcessesList[i_process]));
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteMethod(const std::string& rNameMethod)
{
    for (auto& process : mProcesses)
        process.attr(rNameMethod.c_str())();
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteInitialize()
{
    ExecuteMethod("ExecuteInitialize");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteBeforeSolutionLoop()
{
    ExecuteMethod("ExecuteBeforeSolutionLoop");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteInitializeSolutionStep()
{
    ExecuteMethod("ExecuteInitializeSolutionStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteFinalizeSolutionStep()
{
    ExecuteMethod("ExecuteFinalizeSolutionStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteBeforeOutputStep()
{
    ExecuteMethod("ExecuteBeforeOutputStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteAfterOutputStep()
{
    ExecuteMethod("ExecuteAfterOutputStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::ExecuteFinalize()
{
    ExecuteMethod("ExecuteFinalize");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::IsOutputStep()
{
    ExecuteMethod("IsOutputStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::PrintOutput()
{
    ExecuteMethod("PrintOutput");
}

/***********************************************************************************/
/***********************************************************************************/

void ProcessFactoryUtility::Clear()
{
    ExecuteMethod("Clear");
}
}  // namespace Kratos.
