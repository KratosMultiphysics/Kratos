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
#include "custom_factories/dense_linear_solver_factory.h"

#include "custom_solvers/eigen_sparse_cg_solver.h"
#include "custom_solvers/eigen_sparse_lu_solver.h"
#include "custom_solvers/eigen_sparse_qr_solver.h"
#include "custom_solvers/eigen_direct_solver.h"

#if defined USE_EIGEN_MKL
#include "custom_solvers/eigen_pardiso_lu_solver.h"
#include "custom_solvers/eigen_pardiso_llt_solver.h"
#include "custom_solvers/eigen_pardiso_ldlt_solver.h"
#endif

namespace Kratos
{

void KratosEigenSolversApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosEigenSolversApplication..." << std::endl;

    RegisterDenseLinearSolvers();

    using complex = std::complex<double>;

    // Sparse LU Solver
    using SparseLUType = EigenDirectSolver<EigenSparseLUSolver<double>>;
    static auto SparseLUFactory = SparseLUType::Factory();
    KRATOS_REGISTER_LINEAR_SOLVER("sparse_lu", SparseLUFactory);

    // Complex Sparse LU Solver
    using ComplexSparseLUType = EigenDirectSolver<EigenSparseLUSolver<complex>>;
    static auto ComplexSparseLUFactory = ComplexSparseLUType::Factory();
    KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("sparse_lu_complex", ComplexSparseLUFactory);

    // Sparse QR Solver
    using SparseQRType = EigenDirectSolver<EigenSparseQRSolver<double>>;
    static auto SparseQRFactory = SparseQRType::Factory();
    KRATOS_REGISTER_LINEAR_SOLVER("sparse_qr", SparseQRFactory);

    // Sparse CG Solver
    using SparseCGType = EigenDirectSolver<EigenSparseCGSolver<double>>;
    static auto SparseCGFactory = SparseCGType::Factory();
    KRATOS_REGISTER_LINEAR_SOLVER("sparse_cg", SparseCGFactory);

#if defined USE_EIGEN_MKL

    // Pardiso LU Solver
    using PardisoLUType = EigenDirectSolver<EigenPardisoLUSolver<double>>;
    static auto PardisoLUFactory = PardisoLUType::Factory();
    KRATOS_REGISTER_LINEAR_SOLVER("pardiso_lu", PardisoLUFactory);

    // Pardiso LDLT Solver
    using PardisoLDLTType = EigenDirectSolver<EigenPardisoLDLTSolver<double>>;
    static auto PardisoLDLTFactory = PardisoLDLTType::Factory();
    KRATOS_REGISTER_LINEAR_SOLVER("pardiso_ldlt", PardisoLDLTFactory);

    // Pardiso LLT Solver
    using PardisoLLTType = EigenDirectSolver<EigenPardisoLLTSolver<double>>;
    static auto PardisoLLTFactory = PardisoLLTType::Factory();
    KRATOS_REGISTER_LINEAR_SOLVER("pardiso_llt", PardisoLLTFactory);

    // Complex Pardiso LU Solver
    using ComplexPardisoLUType = EigenDirectSolver<EigenPardisoLUSolver<complex>>;
    static auto ComplexPardisoLUFactory = ComplexPardisoLUType::Factory();
    KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("pardiso_lu_complex", ComplexPardisoLUFactory);

    // Complex Pardiso LDLT Solver
    using ComplexPardisoLDLTType = EigenDirectSolver<EigenPardisoLDLTSolver<complex>>;
    static auto ComplexPardisoLDLTFactory = ComplexPardisoLDLTType::Factory();
    KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("pardiso_ldlt_complex", ComplexPardisoLDLTFactory);

    // Complex Pardiso LLT Solver
    using ComplexPardisoLLTType = EigenDirectSolver<EigenPardisoLLTSolver<complex>>;
    static auto ComplexPardisoLLTFactory = ComplexPardisoLLTType::Factory();
    KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("pardiso_llt_complex", ComplexPardisoLLTFactory);

#endif // defined USE_EIGEN_MKL
}

} // namespace Kratos
