//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_MPI_MANAGER_H_INCLUDED
#define KRATOS_MPI_MANAGER_H_INCLUDED

#include "includes/parallel_environment.h"

namespace Kratos
{
///@addtogroup KratosMPICore
///@{

///@name Kratos Classes
///@{

/// Helper class to manage the MPI lifecycle.
/** This class initializes MPI on construction and finalizes it
 *  on destruction (with appropriate checks for multiple 
 *  initialization or finalization). This object is instantiated
 *  the first time it is needed (as of now, on the first 
 *  MPIDataCommunicator construction) and held by ParallelEnvironment
 *  until the end of the program run.
 */
class KRATOS_API(KRATOS_MPI_CORE) MPIManager: public EnvironmentManager
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::unique_ptr<MPIManager> Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    MPIManager(MPIManager& rOther) = delete;

    /// Destruct the manager, finalizing MPI in the process.
    ~MPIManager() override;

    ///@}
    ///@name Inquiry
    ///@{

    /// Query MPI initialization status.
    /** returns false if MPI_Initialized would return 0, true otherwise. */
    bool IsInitialized() const override;

    /// Query MPI finalization status.
    /** returns false if MPI_Finalized would return 0, true otherwise. */
    bool IsFinalized() const override;

    ///@}
private:
    ///@name Friends
    ///@{

    friend class MPIDataCommunicator;

    ///@}
    ///@name Private LifeCycle
    ///@{
    
    MPIManager();

    ///@}
    ///@name Private Operations
    ///@{

    /// Create a MPIManager instance.
    /** This initializes MPI if it is not initialized yet. */
    static MPIManager::Pointer Create();

    ///@}
};

///@}

///@}
}

#endif