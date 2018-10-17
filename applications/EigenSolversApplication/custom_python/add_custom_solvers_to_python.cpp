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
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigensystem_solver.h"

#include "includes/standard_linear_solver_factory.h"
#include "eigen_solvers_application.h"

namespace Kratos
{

namespace Python
{

template <typename SolverType>
void register_solver(pybind11::module& m, const std::string& name)
{
    using namespace pybind11;

    using Base = DirectSolver<typename SolverType::TGlobalSpace, typename SolverType::TLocalSpace>;

    using EigenDirectSolverType = EigenDirectSolver<SolverType>;

    class_<EigenDirectSolverType, typename EigenDirectSolverType::Pointer, Base >
        (m, name.c_str())
        .def(init<>())
        .def(init<Parameters>())
    ;
}

template <typename SolverType>
void register_eigensystem_solver(pybind11::module& m, const std::string& name = SolverType::Name)
{
    using namespace pybind11;

    using Base = LinearSolver<typename SolverType::TGlobalSpace, typename SolverType::TLocalSpace>;

    using EigenSystemSolverType = EigensystemSolver<SolverType>;

    class_<EigenSystemSolverType, typename EigenSystemSolverType::Pointer, Base >
        (m, name.c_str())
        .def(init<Parameters>())
    ;
}

void AddCustomSolversToPython(pybind11::module& m)
{
    using namespace pybind11;

    using complex = std::complex<double>;

    // --- direct solvers

    register_solver<SparseLU<double>>(m, "SparseLUSolver");
    register_solver<SparseLU<complex>>(m, "ComplexSparseLUSolver");

    register_solver<SparseQR<double>>(m, "SparseQRSolver");
    register_solver<SparseQR<complex>>(m, "ComplexSparseQRSolver");

    #if defined USE_EIGEN_MKL
    // The commented complex solvers need to be tested first
    register_solver<PardisoLLT<double>>(m, "PardisoLLTSolver");
    // register_solver<PardisoLLT<complex>>(m, "ComplexPardisoLLTSolver");

    register_solver<PardisoLDLT<double>>(m, "PardisoLDLTSolver");
    // register_solver<PardisoLDLT<complex>>(m, "ComplexPardisoLDLTSolver");

    register_solver<PardisoLU<double>>(m, "PardisoLUSolver");
    register_solver<PardisoLU<complex>>(m, "ComplexPardisoLUSolver");
    #endif // defined USE_EIGEN_MKL

    // --- eigensystem solver

    #if !defined USE_EIGEN_MKL
    register_eigensystem_solver<SparseLU<double>>(m, "EigensystemSolver");
    #else  // !defined USE_EIGEN_MKL
    register_eigensystem_solver<PardisoLDLT<double>>(m, "EigensystemSolver");
    #endif // !defined USE_EIGEN_MKL
}

} // namespace Python

//Must put this definition here to avoid a problem with multiply defined symbols when including the external C libraries
EigenSolversApplicationRegisterLinearSolvers::EigenSolversApplicationRegisterLinearSolvers()
{
    using complex = std::complex<double>;

    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;
    typedef TUblasSparseSpace<complex> ComplexSpaceType;
    typedef TUblasDenseSpace<complex> ComplexLocalSpaceType;
    //typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef EigenDirectSolver<SparseLU<double>> SparseLUType;
    typedef EigenDirectSolver<SparseLU<complex>> ComplexSparseLUType;
    typedef EigenDirectSolver<SparseQR<double>> SparseQRType;
    typedef EigenDirectSolver<SparseQR<complex>> ComplexSparseQRType;

    //REGISTERING SOLVERS
    typedef LinearSolverFactory<SpaceType,  LocalSpaceType> LinearSolverFactoryType;
    typedef LinearSolverFactory<ComplexSpaceType,  ComplexLocalSpaceType> ComplexLinearSolverFactoryType;

    static auto SparseLUFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SparseLUType>();
    static auto ComplexSparseLUFactory= StandardLinearSolverFactory<ComplexSpaceType,ComplexLocalSpaceType,ComplexSparseLUType>();
    static auto SparseQRFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SparseQRType>();
    static auto ComplexSparseQRFactory= StandardLinearSolverFactory<ComplexSpaceType,ComplexLocalSpaceType,ComplexSparseQRType>();

    KratosComponents<LinearSolverFactoryType>::Add("SparseLUSolver", SparseLUFactory);
    KratosComponents<LinearSolverFactoryType>::Add("eigen_sparse_lu", SparseLUFactory);  // NOTE: Retrocompatibility name
    KratosComponents<ComplexLinearSolverFactoryType>::Add("ComplexSparseLUSolver", ComplexSparseLUFactory);
    KratosComponents<ComplexLinearSolverFactoryType>::Add("complex_eigen_sparse_lu", ComplexSparseLUFactory);  // NOTE: Retrocompatibility name
    KratosComponents<LinearSolverFactoryType>::Add("SparseQRSolver", SparseQRFactory);
    KratosComponents<ComplexLinearSolverFactoryType>::Add("ComplexSparseQRSolver", ComplexSparseQRFactory);

#ifdef USE_EIGEN_MKL
    typedef EigenDirectSolver<PardisoLLT<double>> PardisoLLTSolver;
//     typedef EigenDirectSolver<PardisoLLT<complex>> ComplexPardisoLLTSolver;
    typedef EigenDirectSolver<PardisoLDLT<double>> PardisoLDLTType;
//     typedef EigenDirectSolver<PardisoLDLT<complex>> ComplexPardisoLDLTType;
    typedef EigenDirectSolver<PardisoLU<double>> PardisoLUType;
    typedef EigenDirectSolver<PardisoLU<complex>> ComplexPardisoLUType;

    static auto PardisoLLTFactor= StandardLinearSolverFactory<SpaceType,LocalSpaceType,PardisoLLTSolver>();
//     static auto ComplexPardisoLLTFactory= StandardLinearSolverFactory<ComplexSpaceType,ComplexLocalSpaceType,ComplexPardisoLLTSolver>();
    static auto PardisoLDLTFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,PardisoLDLTType>();
//     static auto ComplexPardisoLDLTFactory= StandardLinearSolverFactory<ComplexSpaceType,ComplexLocalSpaceType,ComplexPardisoLDLTType>();
    static auto PardisoLUFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,PardisoLUType>();
    static auto ComplexPardisoLUFactory= StandardLinearSolverFactory<ComplexSpaceType,ComplexLocalSpaceType,ComplexPardisoLUType>();

    KratosComponents<LinearSolverFactoryType>::Add("PardisoLLTSolver", PardisoLLTFactor);
    KratosComponents<LinearSolverFactoryType>::Add("eigen_pardiso_llt", PardisoLLTFactor); // NOTE: Retrocompatibility name
//     KratosComponents<ComplexLinearSolverFactoryType>::Add("ComplexPardisoLLTSolver", ComplexPardisoLLTFactory);
    KratosComponents<LinearSolverFactoryType>::Add("PardisoLDLTSolver", PardisoLDLTFactory);
    KratosComponents<LinearSolverFactoryType>::Add("eigen_pardiso_ldlt", PardisoLDLTFactory); // NOTE: Retrocompatibility name
//     KratosComponents<ComplexLinearSolverFactoryType>::Add("ComplexPardisoLDLTSolver", ComplexPardisoLDLTFactory);
    KratosComponents<LinearSolverFactoryType>::Add("PardisoLUSolver", PardisoLUFactory);
    KratosComponents<LinearSolverFactoryType>::Add("eigen_pardiso_lu", PardisoLUFactory); // NOTE: Retrocompatibility name
    KratosComponents<ComplexLinearSolverFactoryType>::Add("ComplexPardisoLUSolver", ComplexPardisoLUFactory);
    KratosComponents<ComplexLinearSolverFactoryType>::Add("complex_eigen_pardiso_lu", ComplexPardisoLUFactory); // NOTE: Retrocompatibility name
#endif
}

} // namespace Kratos
