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

#ifndef KRATOS_MPI_ENVIRONMENT_H_INCLUDED
#define KRATOS_MPI_ENVIRONMENT_H_INCLUDED

#include "mpi.h"
#include "includes/data_communicator.h"

namespace Kratos {

/// Helper class to manage the MPI lifecycle
class MPIEnvironment {
public:

    static MPIEnvironment& Instance();

    /// Execute MPI initialization operations.
    /** No MPI operations should be called before this point. */
    void Initialize();

    /// Execute MPI finialization operations.
    /** No MPI operations should be called after this point. */
    void Finalize();

    /// Query MPI initialization status.
    /** returns false if MPI_Initialized would return 0, true otherwise. */
    bool IsInitialized();

    /// Query MPI finalization status.
    /** returns false if MPI_Finalized would return 0, true otherwise. */
    bool IsFinalized();

private:

    MPIEnvironment();

    ~MPIEnvironment();

    MPIEnvironment(const MPIEnvironment&) = delete;

    MPIEnvironment& operator=(const MPIEnvironment) = delete;

    static void Create();

    static MPIEnvironment* mpInstance;

    static bool mDestroyed;

};

}

#endif // KRATOS_MPI_ENVIRONMENT_H_INCLUDED