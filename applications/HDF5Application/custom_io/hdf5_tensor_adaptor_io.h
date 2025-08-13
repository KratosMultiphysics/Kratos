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
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "expression/expression.h"
#include "expression/container_expression.h"

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
class KRATOS_API(HDF5_APPLICATION) TensorAdaptorIO
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TensorAdaptorIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TensorAdaptorIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void Write(
        const std::string& rExpressionName,
        const Expression& rExpression,
        const Parameters Attributes);

    std::pair<Expression::Pointer, Parameters> Read(const std::string& rExpressionName);

    template<class TContainerType>
    void Write(
        const std::string& rContainerExpressionName,
        const ContainerExpression<TContainerType>& rContainerExpression,
        const Parameters Attributes);

    template<class TContainerType>
    Parameters Read(
        const std::string& rExpressionName,
        ContainerExpression<TContainerType>& rContainerExpression);

    ///@}

private:
    ///@name Private member variables
    ///@{

    File::Pointer mpFile;

    std::string mPrefix;

    ///@}

}; // class TensorAdaptorIO.


///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.
