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

#if !defined(KRATOS_HDF5_NODAL_SOLUTION_STEP_BOSSAK_IO_H_INCLUDED)
#define KRATOS_HDF5_NODAL_SOLUTION_STEP_BOSSAK_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_io/hdf5_nodal_solution_step_data_io.h"

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

/// A class for IO of nodal solution step data in HDF5 with weighted Bossak acceleration.
/**
 * This class performs regular IO of nodal solution step data except for the
 * ACCELERATION, which is stored as a weighted combination of the current and
 * previous time step according to the Bossak scheme.
 * 
 */
class NodalSolutionStepBossakIO : private NodalSolutionStepDataIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NodalSolutionStepBossakIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalSolutionStepBossakIO(Parameters Settings, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void WriteNodalResults(NodesContainerType const& rNodes);

    void ReadNodalResults(NodesContainerType& rNodes, Communicator& rComm);

    void SetAlphaBossak(double alpha) noexcept;

    ///@}

private:
    ///@name Member Variables
    ///@{
    double mAlphaBossak = -0.3;
    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_SOLUTION_STEP_BOSSAK_IO_H_INCLUDED defined
