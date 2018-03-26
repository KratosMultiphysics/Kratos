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

#if !defined(KRATOS_HDF5_NON_HISTORICAL_NODAL_VALUE_IO_H_INCLUDED)
#define KRATOS_HDF5_NON_HISTORICAL_NODAL_VALUE_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

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
class NonHistoricalNodalValueIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NonHistoricalNodalValueIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NonHistoricalNodalValueIO(Parameters Settings, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void WriteNodalResults(NodesContainerType const& rNodes);

    void ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm);
    
    ///@}

private:
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;
    std::string mPrefix;
    std::vector<std::string> mVariableNames;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // class NonHistoricalNodalValueIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NON_HISTORICAL_NODAL_VALUE_IO_H_INCLUDED defined
