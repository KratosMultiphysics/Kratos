//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
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
}

}  // namespace Kratos.


