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

namespace Kratos {

/// Helper utilities to manage the MPI lifecycle
namespace MPIEnvironment {

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

}
}