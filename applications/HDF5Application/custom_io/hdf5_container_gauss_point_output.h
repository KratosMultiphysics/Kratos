//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_HDF5_CONTAINER_GAUSS_POINT_OUTPUT_H_INCLUDED)
#define KRATOS_HDF5_CONTAINER_GAUSS_POINT_OUTPUT_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"
#include "includes/process_info.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "hdf5_application_define.h"

namespace Kratos
{
class Parameters;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of element data in HDF5.
template <typename TContainerType, typename... TComponents>
class ContainerGaussPointOutput
{
public:
    ///@name Type Definitions
    ///@{

    using TContainerItemType = typename TContainerType::value_type;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ContainerGaussPointOutput);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ContainerGaussPointOutput(Parameters Settings, File::Pointer pFile, const std::string& rPath);

    ///@}
    ///@name Operations
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void WriteContainerGaussPointsValues(
        TContainerType& rContainerItems,
        const DataCommunicator& rDataCommunicator,
        const ProcessInfo& rProcessInfo);

    ///@}

private:
    ///@name Member Variables
    ///@{
    File::Pointer mpFile;
    std::string mVariablePath;
    std::vector<std::string> mVariableNames;
    ///@}
    ///@name Private Operations
    ///@{

    template <typename... TArgs>
    void WriteRegisteredGaussPointValues(
        const std::string& rComponentName,
        TArgs&... args);

    ///@}

}; // class ContainerGaussPointOutput.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONTAINER_GAUSS_POINT_OUTPUT_H_INCLUDED defined
