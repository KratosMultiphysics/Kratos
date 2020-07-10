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

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

// Linear solvers
#include "custom_factories/dense_linear_solver_factory.h"
#include "custom_solvers/eigen_dense_direct_solver.h"
#include "custom_solvers/eigen_dense_col_piv_householder_qr_solver.h"
#include "custom_solvers/eigen_dense_householder_qr_solver.h"
#include "custom_solvers/eigen_dense_llt_solver.h"
#include "custom_solvers/eigen_dense_partial_piv_lu_solver.h"

namespace Kratos {

void RegisterDenseLinearSolvers()
{
    using complex = std::complex<double>;

    // Dense ColPivHouseholderQR solver
    using EigenDirectColPivHouseholderQRType = EigenDenseDirectSolver<EigenDenseColPivHouseholderQRSolver<double>>;
    static auto DenseColPivHouseholderQRSolverFactory = EigenDirectColPivHouseholderQRType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_col_piv_householder_qr", DenseColPivHouseholderQRSolverFactory);

    // Dense HouseholderQR solver
    using EigenDirectHouseholderQRType = EigenDenseDirectSolver<EigenDenseHouseholderQRSolver<double>>;
    static auto DenseHouseholderQRSolverFactory = EigenDirectHouseholderQRType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_householder_qr", DenseHouseholderQRSolverFactory);

    // Dense LLT solver
    using EigenDirectLLTType = EigenDenseDirectSolver<EigenDenseLLTSolver<double>>;
    static auto DenseLLTSolverFactory = EigenDirectLLTType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_llt", DenseLLTSolverFactory);

    // Dense PartialPivLU solver
    using EigenDirectPartialPivLUType = EigenDenseDirectSolver<EigenDensePartialPivLUSolver<double>>;
    static auto DensePartialPivLUSolverFactory = EigenDirectPartialPivLUType::Factory();
    KRATOS_REGISTER_DENSE_LINEAR_SOLVER("dense_partial_piv_lu", DensePartialPivLUSolverFactory);

    // Complex dense ColPivHouseholderQR solver
    using ComplexEigenDirectColPivHouseholderQRType = EigenDenseDirectSolver<EigenDenseColPivHouseholderQRSolver<complex>>;
    static auto ComplexDenseColPivHouseholderQRSolverFactory = ComplexEigenDirectColPivHouseholderQRType::Factory();
    KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER("complex_dense_col_piv_householder_qr", ComplexDenseColPivHouseholderQRSolverFactory);

    // Complex dense HouseholderQR solver
    using ComplexEigenDirectHouseholderQRType = EigenDenseDirectSolver<EigenDenseHouseholderQRSolver<complex>>;
    static auto ComplexDenseHouseholderQRSolverFactory = ComplexEigenDirectHouseholderQRType::Factory();
    KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER("complex_dense_householder_qr", ComplexDenseHouseholderQRSolverFactory);

    // Complex dense PartialPivLU solver
    using ComplexEigenDirectPartialPivLUType = EigenDenseDirectSolver<EigenDensePartialPivLUSolver<complex>>;
    static auto ComplexDensePartialPivLUSolverFactory = ComplexEigenDirectPartialPivLUType::Factory();
    KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER("complex_dense_partial_piv_lu", ComplexDensePartialPivLUSolverFactory);

}

template class KratosComponents<DenseLinearSolverFactoryType>;
template class KratosComponents<ComplexDenseLinearSolverFactoryType>;
}