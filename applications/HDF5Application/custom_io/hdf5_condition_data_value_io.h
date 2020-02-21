//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_HDF5_CONDITION_DATA_VALUE_IO_H_INCLUDED)
#define KRATOS_HDF5_CONDITION_DATA_VALUE_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_io/hdf5_container_component_io.h"
#include "custom_io/hdf5_file.h"
#include "hdf5_application_define.h"
#include "includes/communicator.h"

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
class ConditionDataValueIO : public ContainerComponentIO<ConditionsContainerType,
                                                       ConditionType,
                                                       Variable<array_1d<double, 3>>,
                                                       Variable<double>,
                                                       Variable<int>,
                                                       Variable<Vector<double>>,
                                                       Variable<Matrix<double>>>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ContainerComponentIO<ConditionsContainerType,
                                          ConditionType,
                                          Variable<array_1d<double, 3>>,
                                          Variable<double>,
                                          Variable<int>,
                                          Variable<Vector<double>>,
                                          Variable<Matrix<double>>>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConditionDataValueIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ConditionDataValueIO(Parameters Settings, File::Pointer pFile)
        : BaseType(Settings, pFile, "/ConditionDataValues")
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void WriteConditionResults(ConditionsContainerType const& rConditions)
    {
        this->WriteContainerComponents(rConditions);
    }

    void ReadConditionResults(ConditionsContainerType& rConditions, Communicator& rComm)
    {
        this->ReadContainerComponents(rConditions, rComm);
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    ///@}

private:
    ///@name Member Variables
    ///@{
    ///@}
    ///@name Private Operations
    ///@{
    ///@}

}; // class ConditionDataValueIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONDITION_DATA_VALUE_IO_H_INCLUDED defined
