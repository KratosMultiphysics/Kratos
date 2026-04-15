//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "trilinos_application.h"
#include "custom_factories/trilinos_linear_solver_factory.h"

namespace Kratos
{
void KratosTrilinosApplication::Register()
{
    KRATOS_INFO("") << "    KRATOS  _____     _ _ _\n"
                    << "           |_   _| __(_) (_)_ __   ___  ___\n"
                    << "             | || '__| | | | '_ \\ / _ \\/ __|\n"
                    << "             | || |  | | | | | | | (_) \\__ \\\n"
                    << "             |_||_|  |_|_|_|_| |_|\\___/|___/\n"
                    << "Initializing KratosTrilinosApplication..." << std::endl;

    RegisterTrilinosLinearSolvers();

#ifdef HAVE_TPETRA
    // Initialize Kokkos (and Tpetra) via the ScopeGuard so that Tpetra objects
    // can be used for the lifetime of the application.  We pass dummy argc/argv
    // because command-line argument parsing is not relevant here. 
    // TODO: Think a way to allow users to pass command-line arguments to Tpetra/Kokkos if they want to (e.g. for Kokkos backend selection or Tpetra debug options). Maybe environment variables are sufficient for this purpose, and would be simpler than trying to pass command-line arguments through Kratos.
    if (!mpTpetraScope) {
        // int argc = 0;
        // char** argv = nullptr;
        // mpTpetraScope = std::make_unique<Tpetra::ScopeGuard>(&argc, &argv);
    }
#endif
}

void KratosTrilinosApplication::DeregisterApplication()
{
#ifdef HAVE_TPETRA
    // Reset the ScopeGuard, which calls Tpetra::finalize() / Kokkos::finalize()
    // and marks Tpetra as no longer initialised.
    if (mpTpetraScope) {
        mpTpetraScope.reset();
    }
#endif
}

}  // namespace Kratos.


