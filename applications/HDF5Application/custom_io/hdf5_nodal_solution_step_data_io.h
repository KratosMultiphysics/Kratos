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
//

#if !defined(KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED)
#define KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED

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

/// A class for IO of nodal solution step data in HDF5.
class NodalSolutionStepDataIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NodalSolutionStepDataIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalSolutionStepDataIO(Parameters Settings, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void WriteNodalResults(NodesContainerType const& rNodes, unsigned Step=0);

    void ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm, unsigned Step=0);
    
    ///@}

protected:
    ///@name Protected Operations
    ///@{

    std::string const& GetPrefix() const noexcept
    {
        return mPrefix;
    }

    std::vector<std::string> const& VariableNames() const noexcept
    {
        return mVariableNames;
    }

    File& GetFile()
    {
        return *mpFile;
    }

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

}; // class NodalSolutionStepDataIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_SOLUTION_STEP_DATA_IO_H_INCLUDED defined
