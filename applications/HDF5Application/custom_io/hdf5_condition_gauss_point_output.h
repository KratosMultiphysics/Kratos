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

#if !defined(KRATOS_HDF5_CONDITION_GAUSS_POINT_OUTPUT_H_INCLUDED)
#define KRATOS_HDF5_CONDITION_GAUSS_POINT_OUTPUT_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/process_info.h"

// Application includes
#include "custom_io/hdf5_container_gauss_point_output.h"
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
class ConditionGaussPointOutput : public ContainerGaussPointOutput<ConditionsContainerType,
                                                                   Variable<array_1d<double, 3>>,
                                                                   Variable<double>,
                                                                   Variable<int>,
                                                                   Variable<Vector<double>>,
                                                                   Variable<Matrix<double>>>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ContainerGaussPointOutput<ConditionsContainerType,
                                               Variable<array_1d<double, 3>>,
                                               Variable<double>,
                                               Variable<int>,
                                               Variable<Vector<double>>,
                                               Variable<Matrix<double>>>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConditionGaussPointOutput);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ConditionGaussPointOutput(Parameters Settings, File::Pointer pFile)
        : BaseType(Settings, pFile, "/ConditionGaussPointValues")
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void WriteConditionGaussPointValues(
        ConditionsContainerType& rConditions,
        const DataCommunicator& rDataCommunicator,
        const ProcessInfo& rProcessInfo)
    {
        this->WriteContainerGaussPointsValues(rConditions, rDataCommunicator, rProcessInfo);
    }

    ///@}

}; // class ConditionGaussPointOutput.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONDITION_GAUSS_POINT_OUTPUT_H_INCLUDED defined
