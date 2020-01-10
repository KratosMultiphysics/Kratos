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
#include "custom_python/add_custom_solvers_to_python.h"
#include "linear_solvers/linear_solver.h"
#include "custom_solvers/eigen_sparse_cg_solver.h"
#include "custom_solvers/eigen_sparse_lu_solver.h"
#include "custom_solvers/eigen_sparse_qr_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigensystem_solver.h"

#if defined USE_EIGEN_MKL
#include "custom_solvers/eigen_pardiso_lu_solver.h"
#include "custom_solvers/eigen_pardiso_llt_solver.h"
#include "custom_solvers/eigen_pardiso_ldlt_solver.h"
#endif

#include "factories/standard_linear_solver_factory.h"
#include "eigen_solvers_application.h"

namespace Kratos {
namespace Python {

template <typename SolverType>
void register_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using GlobalSpace = typename SpaceType<typename SolverType::Scalar>::Global;
    using LocalSpace = typename SpaceType<typename SolverType::Scalar>::Local;

    using Type = EigenDirectSolver<SolverType>;
    using Holder = typename Type::Pointer;
    using Base = DirectSolver<GlobalSpace, LocalSpace>;

    py::class_<Type, Holder, Base>
        (m, name.c_str())
        .def(py::init<>())
        .def(py::init<Parameters>())
    ;
}

void register_eigensystem_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using Base = LinearSolver<UblasSpace<double, CompressedMatrix, Vector>,
        UblasSpace<double, Matrix, Vector>>;

    using EigenSystemSolverType = EigensystemSolver<>;

    py::class_<EigenSystemSolverType, typename EigenSystemSolverType::Pointer, Base >
        (m, name.c_str())
        .def(py::init<Parameters>())
    ;
}

void AddCustomSolversToPython(pybind11::module& m)
{
    using complex = std::complex<double>;

    // --- direct solvers

    register_solver<EigenSparseLUSolver<double>>(m, "SparseLUSolver");
    register_solver<EigenSparseCGSolver<double>>(m, "SparseCGSolver");
    register_solver<EigenSparseQRSolver<double>>(m, "SparseQRSolver");

    register_solver<EigenSparseLUSolver<complex>>(m, "ComplexSparseLUSolver");

#if defined USE_EIGEN_MKL
    register_solver<EigenPardisoLUSolver<double>>(m, "PardisoLUSolver");
    register_solver<EigenPardisoLDLTSolver<double>>(m, "PardisoLDLTSolver");
    register_solver<EigenPardisoLLTSolver<double>>(m, "PardisoLLTSolver");

    register_solver<EigenPardisoLUSolver<complex>>(m, "ComplexPardisoLUSolver");
    register_solver<EigenPardisoLDLTSolver<complex>>(m, "ComplexPardisoLDLTSolver");
    register_solver<EigenPardisoLLTSolver<complex>>(m, "ComplexPardisoLLTSolver");
#endif // defined USE_EIGEN_MKL

    // --- eigensystem solver

    register_eigensystem_solver(m, "EigensystemSolver");
}

} // namespace Python

//Must put this definition here to avoid a problem with multiply defined symbols when including the external C libraries
EigenSolversApplicationRegisterLinearSolvers::EigenSolversApplicationRegisterLinearSolvers()
{
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
