/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "eigen_solvers_application.h"

namespace Kratos
{
KratosEigenSolversApplication::KratosEigenSolversApplication() : KratosApplication("EigenSolversApplication") {}

void KratosEigenSolversApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosEigenSolversApplication... " << std::endl;

    EigenSolversApplicationRegisterLinearSolvers();
}

} // namespace Kratos
