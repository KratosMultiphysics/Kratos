/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
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
#include "custom_solvers/eigen_dense_eigenvalue_solver.h"
#include "custom_solvers/eigensystem_solver.h"

#if defined USE_EIGEN_MKL
#include "custom_solvers/eigen_pardiso_lu_solver.h"
#include "custom_solvers/eigen_pardiso_llt_solver.h"
#include "custom_solvers/eigen_pardiso_ldlt_solver.h"
#endif

#if defined USE_EIGEN_FEAST
#include "custom_solvers/feast_eigensystem_solver.h"
#endif

#include "custom_solvers/spectra_sym_g_eigs_shift_solver.h"
#include "custom_solvers/spectra_g_eigs_shift_solver.h"


#include "factories/standard_linear_solver_factory.h"

/* Utilities */
#include "custom_utilities/feast_condition_number_utility.h"

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

void register_dense_eigenvalue_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using LocalSpace = typename SpaceType<double>::Local;

    using Type = DenseEigenvalueSolver<>;
    using Holder = typename Type::Pointer;
    using Base = LinearSolver<LocalSpace, LocalSpace>;

    void (Base::*pointer_to_solve_dense)(Base::SparseMatrixType& rA, Base::SparseMatrixType& rDummy, Base::DenseVectorType& rX, Base::DenseMatrixType& rB) = &Base::Solve;

    py::class_<Type, Holder, Base>
        (m, name.c_str())
        .def(py::init<Parameters>())
        .def("Solve",pointer_to_solve_dense)
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

template<typename EigenSystemSolverType>
void register_feast_eigensystem_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using DataTypeIn = typename EigenSystemSolverType::ValueTypeIn;
    using DataTypeOut = typename EigenSystemSolverType::ValueTypeOut;
    using SparseSpaceType = TUblasSparseSpace<DataTypeIn>;
    using DenseSpaceType = TUblasDenseSpace<DataTypeOut>;
    using Base = LinearSolver<SparseSpaceType, DenseSpaceType>;

    py::class_<EigenSystemSolverType, typename EigenSystemSolverType::Pointer, Base >
        (m, name.c_str())
        .def(py::init<Parameters>())
    ;
}

void register_spectra_g_eigs_shift_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using Base = LinearSolver<UblasSpace<double, CompressedMatrix, Vector>,
        UblasSpace<double, Matrix, Vector>>;

    using SpectraGEigsRealSolverType = SpectraGEigsShiftSolver<>;

    py::class_<SpectraGEigsRealSolverType, typename SpectraGEigsRealSolverType::Pointer, Base >
        (m, name.c_str())
        .def(py::init<Parameters>())
    ;
}

void register_spectra_sym_g_eigs_shift_solver(pybind11::module& m, const std::string& name)
{
    namespace py = pybind11;

    using Base = LinearSolver<UblasSpace<double, CompressedMatrix, Vector>,
        UblasSpace<double, Matrix, Vector>>;

    using SpectraSymGEigsRealSolverType = SpectraSymGEigsShiftSolver<>;

    py::class_<SpectraSymGEigsRealSolverType, typename SpectraSymGEigsRealSolverType::Pointer, Base >
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
    bool (DenseLinearSolverType::*pointer_to_multi_solve_dense)(DenseLinearSolverType::SparseMatrixType& rA, DenseLinearSolverType::DenseMatrixType& rX, DenseLinearSolverType::DenseMatrixType& rB) = &DenseLinearSolverType::Solve;
    bool (ComplexDenseLinearSolverType::*pointer_to_complex_solve_dense)(ComplexDenseLinearSolverType::SparseMatrixType& rA, ComplexDenseLinearSolverType::VectorType& rX, ComplexDenseLinearSolverType::VectorType& rB) = &ComplexDenseLinearSolverType::Solve;

    py::class_<DenseLinearSolverType, DenseLinearSolverType::Pointer>(m,"DenseLinearSolver")
        .def(py::init< >())
        .def("Initialize",&DenseLinearSolverType::Initialize)
        .def("Solve",pointer_to_solve_dense)
        .def("Solve",pointer_to_multi_solve_dense)
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
    namespace py = pybind11;

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

    // --- dense eigenvalue solver
    register_dense_eigenvalue_solver(m, "DenseEigenvalueSolver");

#if defined USE_EIGEN_FEAST
    register_feast_eigensystem_solver<FEASTEigensystemSolver<true, double, double>>(m, "FEASTSymmetricEigensystemSolver");
    register_feast_eigensystem_solver<FEASTEigensystemSolver<false, double, complex>>(m, "FEASTGeneralEigensystemSolver");
    register_feast_eigensystem_solver<FEASTEigensystemSolver<true, complex, complex>>(m, "ComplexFEASTSymmetricEigensystemSolver");
    register_feast_eigensystem_solver<FEASTEigensystemSolver<false, complex, complex>>(m, "ComplexFEASTGeneralEigensystemSolver");
#endif

    // --- spectra eigensystem solver
    register_spectra_sym_g_eigs_shift_solver(m, "SpectraSymGEigsShiftSolver") ;
    register_spectra_g_eigs_shift_solver(m, "SpectraGEigsShiftSolver");

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef FEASTConditionNumberUtility<SparseSpaceType, LocalSpaceType> FEASTConditionNumberUtilityType;
    py::class_<FEASTConditionNumberUtilityType,FEASTConditionNumberUtilityType::Pointer>(m,"FEASTConditionNumberUtility")
        .def("GetConditionNumber", &FEASTConditionNumberUtilityType::GetConditionNumber)
        ;
}

} // namespace Python

} // namespace Kratos
