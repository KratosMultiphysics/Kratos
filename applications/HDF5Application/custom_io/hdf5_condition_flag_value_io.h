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

#if !defined(KRATOS_HDF5_CONDITION_FLAG_VALUE_IO_H_INCLUDED)
#define KRATOS_HDF5_CONDITION_FLAG_VALUE_IO_H_INCLUDED

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
class ConditionFlagValueIO
    : public ContainerComponentIO<ConditionsContainerType, ConditionType, Flags>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ContainerComponentIO<ConditionsContainerType, ConditionType, Flags>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConditionFlagValueIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ConditionFlagValueIO(Parameters Settings, File::Pointer pFile)
        : BaseType(Settings, pFile, "/ConditionFlagValues")
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void WriteConditionFlags(ConditionsContainerType const& rConditions)
    {
        this->WriteContainerComponents(rConditions);
    }

    void ReadConditionFlags(ConditionsContainerType& rConditions, Communicator& rComm)
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

}; // class ConditionFlagValueIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONDITION_FLAG_VALUE_IO_H_INCLUDED defined
