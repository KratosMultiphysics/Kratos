//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#pragma once

// System includes
#include <string>
#include <tuple>
#include <variant>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "containers/nd_data.h"

// Application includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of element data in HDF5.
class KRATOS_API(HDF5_APPLICATION) NDDataIO
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NDDataPointerType = std::variant<
                                            NDData<unsigned char>::Pointer,
                                            NDData<bool>::Pointer,
                                            NDData<int>::Pointer,
                                            NDData<double>::Pointer
                                        >;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NDDataIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NDDataIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void Write(
        const std::string& rNDDataName,
        NDDataPointerType pNDData,
        const Parameters Attributes);

    std::pair<NDDataPointerType, Parameters> Read(const std::string& rNDDataName);

    ///@}

private:
    ///@name Private member variables
    ///@{

    File::Pointer mpFile;

    std::string mPrefix;

    ///@}

}; // class NDDataIO.


///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.
