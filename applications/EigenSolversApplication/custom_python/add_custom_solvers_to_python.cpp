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
#include "includes/define_python.h"
#include "custom_python/add_custom_solvers_to_python.h"
#include "linear_solvers/linear_solver.h"
#include "custom_solvers/eigen_sparse_cg_solver.h"
#include "custom_solvers/eigen_sparse_lu_solver.h"
#include "custom_solvers/eigen_sparse_qr_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigen_dense_col_piv_householder_qr_solver.h"
#include "custom_solvers/eigen_dense_householder_qr_solver.h"
#include "custom_solvers/eigen_dense_llt_solver.h"
#include "custom_solvers/eigen_dense_partial_piv_lu_solver.h"
#include "custom_solvers/eigen_dense_direct_solver.h"
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

template<typename SolverType>
void register_dense_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using LocalSpace = typename SpaceType<typename SolverType::Scalar>::Local;

    using Type = EigenDenseDirectSolver<SolverType>;
    using Holder = typename Type::Pointer;
    using Base = DirectSolver<LocalSpace, LocalSpace>;

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

void register_base_dense_solver(pybind11::module& m)
{
    namespace py = pybind11;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;
    typedef LinearSolver<LocalSpaceType,  LocalSpaceType> DenseLinearSolverType;
    typedef LinearSolver<ComplexLocalSpaceType, ComplexLocalSpaceType> ComplexDenseLinearSolverType;

    bool (DenseLinearSolverType::*pointer_to_solve_dense)(DenseLinearSolverType::SparseMatrixType& rA, DenseLinearSolverType::VectorType& rX, DenseLinearSolverType::VectorType& rB) = &DenseLinearSolverType::Solve;
    bool (ComplexDenseLinearSolverType::*pointer_to_complex_solve_dense)(ComplexDenseLinearSolverType::SparseMatrixType& rA, ComplexDenseLinearSolverType::VectorType& rX, ComplexDenseLinearSolverType::VectorType& rB) = &ComplexDenseLinearSolverType::Solve;

    py::class_<DenseLinearSolverType, DenseLinearSolverType::Pointer>(m,"DenseLinearSolver")
        .def(py::init< >())
        .def("Initialize",&DenseLinearSolverType::Initialize)
        .def("Solve",pointer_to_solve_dense)
        .def("Clear",&DenseLinearSolverType::Clear)
        .def("__str__",PrintObject<DenseLinearSolverType>)
        ;

    py::class_<ComplexDenseLinearSolverType, ComplexDenseLinearSolverType::Pointer>(m,"ComplexDenseLinearSolver")
        .def(py::init< >())
        .def("Initialize",&ComplexDenseLinearSolverType::Initialize)
        .def("Solve",pointer_to_complex_solve_dense)
        .def("Clear",&ComplexDenseLinearSolverType::Clear)
        .def("__str__",PrintObject<ComplexDenseLinearSolverType>)
        ;

    typedef DirectSolver<LocalSpaceType,  LocalSpaceType> DenseDirectSolverType;
    py::class_<DenseDirectSolverType, DenseDirectSolverType::Pointer, DenseLinearSolverType>(m,"DirectSolver")
        .def(py::init< >() )
        .def(py::init<Parameters>())
        .def("__str__", PrintObject<DenseDirectSolverType>)
        ;

    typedef DirectSolver<ComplexLocalSpaceType,  ComplexLocalSpaceType> ComplexDenseDirectSolverType;
    py::class_<ComplexDenseDirectSolverType, ComplexDenseDirectSolverType::Pointer, ComplexDenseLinearSolverType>(m,"ComplexDirectSolver")
    .def(py::init< >() )
    .def(py::init<Parameters>())
    .def("__str__", PrintObject<ComplexDenseDirectSolverType>)
    ;

    typedef LinearSolverFactory< LocalSpaceType, LocalSpaceType > DenseLinearSolverFactoryType;
    py::class_<DenseLinearSolverFactoryType, DenseLinearSolverFactoryType::Pointer>(m, "DenseLinearSolverFactory")
        .def( py::init< >() )
        .def("Create",&DenseLinearSolverFactoryType::Create)
        .def("Has",&DenseLinearSolverFactoryType::Has)
        ;

    typedef LinearSolverFactory< ComplexLocalSpaceType, ComplexLocalSpaceType > ComplexDenseLinearSolverFactoryType;
    py::class_<ComplexDenseLinearSolverFactoryType, ComplexDenseLinearSolverFactoryType::Pointer>(m, "ComplexDenseLinearSolverFactory")
        .def( py::init< >() )
        .def("Create",&ComplexDenseLinearSolverFactoryType::Create)
        .def("Has",&ComplexDenseLinearSolverFactoryType::Has)
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

    register_base_dense_solver(m);

    register_dense_solver<EigenDenseColPivHouseholderQRSolver<double>>(m, "DenseColPivHouseholderQRSolver");
    register_dense_solver<EigenDenseHouseholderQRSolver<double>>(m, "DenseHouseholderQRSolver");
    register_dense_solver<EigenDenseLLTSolver<double>>(m, "DenseLLTSolver");
    register_dense_solver<EigenDensePartialPivLUSolver<double>>(m, "DensePartialPivLUSolver");

    register_dense_solver<EigenDenseColPivHouseholderQRSolver<complex>>(m, "ComplexDenseColPivHouseholderQRSolver");
    register_dense_solver<EigenDenseHouseholderQRSolver<complex>>(m, "ComplexDenseHouseholderQRSolver");
    register_dense_solver<EigenDensePartialPivLUSolver<complex>>(m, "ComplexDensePartialPivLUSolver");

    // --- eigensystem solver

    register_eigensystem_solver(m, "EigensystemSolver");
}

} // namespace Python

} // namespace Kratos
