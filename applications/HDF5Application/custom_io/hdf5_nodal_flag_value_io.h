//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    https://github.com/msandre
//

#if !defined(KRATOS_HDF5_NODAL_FLAG_VALUE_IO_H_INCLUDED)
#define KRATOS_HDF5_NODAL_FLAG_VALUE_IO_H_INCLUDED

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

namespace Kratos
{
class Parameters;
class Communicator;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of non-historical nodal values in HDF5.
class NodalFlagValueIO : public ContainerComponentIO<NodesContainerType, NodeType, Flags>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ContainerComponentIO<NodesContainerType, NodeType, Flags>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NodalFlagValueIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalFlagValueIO(Parameters Settings, File::Pointer pFile)
        : BaseType(Settings, pFile, "/NodalFlagValues")
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void WriteNodalFlags(NodesContainerType const& rNodes)
    {
        this->WriteContainerComponents(rNodes);
    }

    void ReadNodalFlags(NodesContainerType& rNodes, Communicator& rComm)
    {
        this->ReadContainerComponents(rNodes, rComm);
        rComm.SynchronizeNodalFlags();
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // class NodalFlagValueIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_FLAG_VALUE_IO_H_INCLUDED defined
