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

#if !defined(KRATOS_HDF5_NODAL_SOLUTION_STEP_VARIABLES_IO_H_INCLUDED)
#define KRATOS_HDF5_NODAL_SOLUTION_STEP_VARIABLES_IO_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/variables_list.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class for IO of nodal solution step variables in HDF5.
class NodalSolutionStepVariablesIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NodalSolutionStepVariablesIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalSolutionStepVariablesIO(std::string Prefix, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void WriteVariablesList(ModelPart const& rModelPart);

    void ReadAndAssignVariablesList(ModelPart& rModelPart) const;

    void WriteBufferSize(int BufferSize);

    void ReadAndAssignBufferSize(ModelPart& rModelPart) const;

    ///@}

private:
    ///@name Member Variables
    ///@{
    std::string mPrefix;
    File::Pointer mpFile;
    ///@}

}; // class NodalSolutionStepVariablesIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_SOLUTION_STEP_VARIABLES_IO_H_INCLUDED defined
