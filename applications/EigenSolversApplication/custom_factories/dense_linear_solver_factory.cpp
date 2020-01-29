/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Quirin Aumann
*/

// // Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

// // Linear solvers
#include "custom_solvers/eigen_dense_direct_solver.h"
#include "custom_solvers/eigen_dense_colpivhouseholderqr_solver.h"
#include "custom_solvers/eigen_dense_householderqr_solver.h"
#include "custom_solvers/eigen_dense_llt_solver.h"
#include "custom_solvers/eigen_dense_partialpivlu_solver.h"

namespace Kratos {

void RegisterDenseLinearSolvers()
{
    using complex = std::complex<double>;

    // Dense ColPivHouseholderQR solver
    using EigenDirectColPivHouseholderQRType = EigenDenseDirectSolver<EigenDenseColPivHouseholderQRSolver<double>>;
    static auto DenseColPivHouseholderQRSolverFactory = EigenDirectColPivHouseholderQRType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_colpivhouseholderqr", DenseColPivHouseholderQRSolverFactory);

    // Dense HouseholderQR solver
    using EigenDirectHouseholderQRType = EigenDenseDirectSolver<EigenDenseHouseholderQRSolver<double>>;
    static auto DenseHouseholderQRSolverFactory = EigenDirectHouseholderQRType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_householderqr", DenseHouseholderQRSolverFactory);

    // Dense LLT solver
    using EigenDirectLLTType = EigenDenseDirectSolver<EigenDenseLLTSolver<double>>;
    static auto DenseLLTSolverFactory = EigenDirectLLTType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_llt", DenseLLTSolverFactory);

    // Dense PartialPivLU solver
    using EigenDirectPartialPivLUType = EigenDenseDirectSolver<EigenDensePartialPivLUSolver<double>>;
    static auto DensePartialPivLUSolverFactory = EigenDirectPartialPivLUType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_partialpivlu", DensePartialPivLUSolverFactory);

    // Complex dense ColPivHouseholderQR solver
    using ComplexEigenDirectColPivHouseholderQRType = EigenDenseDirectSolver<EigenDenseColPivHouseholderQRSolver<complex>>;
    static auto ComplexDenseColPivHouseholderQRSolverFactory = ComplexEigenDirectColPivHouseholderQRType::Factory();
    KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER("complex_dense_colpivhouseholderqr", ComplexDenseColPivHouseholderQRSolverFactory);

    // Complex dense HouseholderQR solver
    using ComplexEigenDirectHouseholderQRType = EigenDenseDirectSolver<EigenDenseHouseholderQRSolver<complex>>;
    static auto ComplexDenseHouseholderQRSolverFactory = ComplexEigenDirectHouseholderQRType::Factory();
    KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER("complex_dense_householderqr", ComplexDenseHouseholderQRSolverFactory);

    // Complex dense PartialPivLU solver
    using ComplexEigenDirectPartialPivLUType = EigenDenseDirectSolver<EigenDensePartialPivLUSolver<complex>>;
    static auto ComplexDensePartialPivLUSolverFactory = ComplexEigenDirectPartialPivLUType::Factory();
    KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER("complex_dense_partialpivlu", ComplexDensePartialPivLUSolverFactory);

}

template class KratosComponents<DenseLinearSolverFactoryType>;
}