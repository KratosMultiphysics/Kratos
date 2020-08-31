//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_TRILINOS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_TRILINOS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "Epetra_MpiComm.h"
#include "mpi.h"

// Application includes
#include "custom_processes/rans_wall_distance_calculation_process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(RANS_APPLICATION) TrilinosRansWallDistanceCalculationProcess
    : public RansWallDistanceCalculationProcess
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = RansWallDistanceCalculationProcess;

    /// Pointer definition of TrilinosRansWallDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosRansWallDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TrilinosRansWallDistanceCalculationProcess(
        Model& rModel,
        Parameters rParameters)
        : BaseType(rModel, rParameters)
    {
    }

    /// Destructor.
    ~TrilinosRansWallDistanceCalculationProcess() override = default;

    /// Assignment operator.
    TrilinosRansWallDistanceCalculationProcess& operator=(TrilinosRansWallDistanceCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    TrilinosRansWallDistanceCalculationProcess(TrilinosRansWallDistanceCalculationProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    std::string Info() const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    Epetra_MpiComm mMPIComm = Epetra_MpiComm(MPI_COMM_WORLD);

    ///@}
    ///@name Private Operations
    ///@{

    Process::Pointer GetWallDistanceCalculationProcess(
        ModelPart& rModelPart,
        Parameters LinearSolverParameters,
        const int MaxIterations) override;


    ///@}

}; // Class TrilinosRansWallDistanceCalculationProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TrilinosRansWallDistanceCalculationProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_TRILINOS_WALL_DISTANCE_CALCULATION_PROCESS_H_INCLUDED defined
