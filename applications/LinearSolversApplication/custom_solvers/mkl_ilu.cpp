//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// External includes
#include "mkl.h" // MKL_INT, dcsrilu0, dcsrilut
#include "mkl_spblas.h" // mkl_sparse_trsv, mkl_sparse_destroy, mkl_sparse_d_create_csr

// Project includes
#include "custom_solvers/mkl_ilu.hpp" // MKLILUSolverBase
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace

// STL includes
#include <vector> // std::vector
#include <array> // std::array
#include <algorithm> // std::fill
#include <memory> // std::unique_ptr


namespace Kratos {


struct MKLCSRAdaptorDestructor {
    void operator()(sparse_matrix_t* pAdaptor)
    {
        if (pAdaptor) {
            mkl_sparse_destroy(*pAdaptor);
            delete pAdaptor;
        }
        pAdaptor = nullptr;
    }
}; // struct MKLCSRAdaptorDestructor


std::string TranslateMKLSparseReturnCode(sparse_status_t Status)
{
    switch (Status) {
        case SPARSE_STATUS_SUCCESS:
            return "";
        case SPARSE_STATUS_NOT_INITIALIZED:
            return "SPARSE_STATUS_NOT_INITIALIZED: The routine encountered an empty handle or matrix array. ";
        case SPARSE_STATUS_ALLOC_FAILED:
            return "SPARSE_STATUS_ALLOC_FAILED: Internal memory allocation failed.";
        case SPARSE_STATUS_INVALID_VALUE:
            return "SPARSE_STATUS_INVALID_VALUE: The input parameters contain an invalid value.";
        case SPARSE_STATUS_EXECUTION_FAILED:
            return "SPARSE_STATUS_EXECUTION_FAILED: Execution failed.";
        case SPARSE_STATUS_INTERNAL_ERROR:
            return "SPARSE_STATUS_INTERNAL_ERROR: An error in algorithm implementation occurred.";
        case SPARSE_STATUS_NOT_SUPPORTED:
            return "SPARSE_STATUS_NOT_SUPPORTED: The requested operation is not supported.";
        default:
            return "unknown";
    } // switch Status
}


template <class TSparse, class TDense>
struct MKLILUSolverBase<TSparse,TDense>::Impl
{
    std::vector<typename TSparse::DataType> bilut;

    std::vector<MKL_INT> ibilut;

    std::vector<MKL_INT> jbilut;

    std::array<MKL_INT,128> ipar;

    std::array<typename TSparse::DataType,128> dpar;

    std::vector<typename TSparse::DataType> mBuffer;

    std::unique_ptr<sparse_matrix_t,MKLCSRAdaptorDestructor> mpLUAdaptor;

    MKL_INT mIterations;

    typename TSparse::DataType mRelaxation;
}; // struct MKLILUSolverBase::Impl


template <class TSparse, class TDense>
MKLILUSolverBase<TSparse,TDense>::MKLILUSolverBase()
    : MKLILUSolverBase(Parameters())
{
}


template <class TSparse, class TDense>
MKLILUSolverBase<TSparse,TDense>::MKLILUSolverBase(Parameters Settings)
    : Base(),
      mpImpl(new Impl)
{
    KRATOS_TRY
    Settings.AddMissingParameters(this->GetDefaultParameters());
    mpImpl->mIterations = Settings["iterations"].Get<int>();
    mpImpl->mRelaxation = Settings["relaxation"].Get<double>();
    KRATOS_CATCH("")

    std::fill(mpImpl->ipar.begin(), mpImpl->ipar.end(), static_cast<MKL_INT>(0));
    std::fill(mpImpl->dpar.begin(), mpImpl->dpar.end(), static_cast<typename TSparse::DataType>(0));

    mpImpl->mpLUAdaptor = nullptr;
}


template <class TSparse, class TDense>
MKLILUSolverBase<TSparse,TDense>::~MKLILUSolverBase() = default;


/// Compute the incomplete LU factorization.
template <class TSparse, class TDense>
void MKLILUSolverBase<TSparse,TDense>::ProvideAdditionalData(typename Base::SparseMatrix& rLhs,
                                                          typename Base::Vector& rSolution,
                                                          typename Base::Vector& rRhs,
                                                          ModelPart::DofsArrayType& rDofSet,
                                                          ModelPart& rModelPart)
{
    KRATOS_TRY
    Base::ProvideAdditionalData(rLhs, rSolution, rRhs, rDofSet, rModelPart);
    KRATOS_CATCH("")

    MKL_INT error = 0;
    KRATOS_TRY
    if (!rLhs.nnz()) return;

    // Set factorization parameters.
    auto& ipar = mpImpl->ipar;
    auto& dpar = mpImpl->dpar;

    ipar[  4]   = 1e2;      //< maximum number of iterations
    ipar[  5]   = 1;        //< output error messages
    ipar[  6]   = 1;        //< output warning messages
    dpar[ 30]   = 1.0;      //< value to replace small diagonal entries with

    [[maybe_unused]] auto [lhs_view, solution_view, rhs_view] = this->MakeSystemView(
        rLhs,
        &*rSolution.begin(),
        (&*rSolution.begin()) + rSolution.size(),
        &*rRhs.begin(),
        (&*rRhs.begin()) + rRhs.size());

    this->Factorize(mpImpl->ibilut,
                    mpImpl->jbilut,
                    mpImpl->bilut,
                    lhs_view,
                    mpImpl->ipar,
                    mpImpl->dpar);

    KRATOS_CATCH("")
    KRATOS_ERROR_IF(error) << "dcsrilut terminated with error code " << error;

    KRATOS_TRY
    mpImpl->mpLUAdaptor.reset(new sparse_matrix_t);
    const auto status = mkl_sparse_d_create_csr(mpImpl->mpLUAdaptor.get(),
                                                SPARSE_INDEX_BASE_ONE,
                                                static_cast<MKL_INT>(rLhs.size1()),
                                                static_cast<MKL_INT>(rLhs.size2()),
                                                mpImpl->ibilut.data(),
                                                mpImpl->ibilut.data() + 1,
                                                mpImpl->jbilut.data(),
                                                mpImpl->bilut.data());

    if (status != SPARSE_STATUS_SUCCESS) {
        KRATOS_ERROR << "mkl_sparse_d_create_csr failed with return code "
                     << TranslateMKLSparseReturnCode(status);
    } // if status != SPARSE_STATUS_SUCCESS
    KRATOS_CATCH("")

    mpImpl->mBuffer.resize(rRhs.size());
}


template <class TSparse, class TDense>
Parameters MKLILUSolverBase<TSparse,TDense>::GetDefaultParameters()
{
    return Parameters(R"({
        "iterations" : 1,
        "relaxation" : 1.0
    })");
}


template <class TSparse, class TDense>
void MKLILUSolverBase<TSparse,TDense>::Clear()
{
    Base::Clear();
    mpImpl->mpLUAdaptor.reset();
    mpImpl->bilut = std::vector<typename TSparse::DataType>();
    mpImpl->ibilut = std::vector<MKL_INT>();
    mpImpl->jbilut = std::vector<MKL_INT>();
    mpImpl->mBuffer = std::vector<typename TSparse::DataType>();
}


/// Use the factorization and a triangular solver.
template <class TSparse, class TDense>
bool MKLILUSolverBase<TSparse,TDense>::Solve(typename Base::CSRView LhsView,
                                          typename Base::template VectorView</*IsMutable=*/true> SolutionView,
                                          typename Base::template VectorView</*IsMutable=*/false> RhsView)
{
    KRATOS_ERROR_IF_NOT(mpImpl->mpLUAdaptor) << "MKLILUSolverBase::ProvideAdditionalData must be called before MKLILUSolverBase::Solve";
    KRATOS_ERROR_IF_NOT(static_cast<MKL_INT>(mpImpl->mBuffer.size()) == LhsView.row_count);
    KRATOS_ERROR_IF_NOT(LhsView.row_count == LhsView.column_count);
    KRATOS_ERROR_IF_NOT(LhsView.row_count == SolutionView.size);
    KRATOS_ERROR_IF_NOT(SolutionView.size == RhsView.size);

    matrix_descr lu_properties;
    lu_properties.type = SPARSE_MATRIX_TYPE_TRIANGULAR;

    for (MKL_INT i_sweep=0; i_sweep<mpImpl->mIterations; ++i_sweep) {
        // Forward substitution.
        KRATOS_TRY
        lu_properties.mode = SPARSE_FILL_MODE_LOWER;
        lu_properties.diag = SPARSE_DIAG_UNIT;
        const auto status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE,
                                              mpImpl->mRelaxation,
                                              *mpImpl->mpLUAdaptor,
                                              lu_properties,
                                              RhsView.it_begin,
                                              mpImpl->mBuffer.data());
        KRATOS_ERROR_IF_NOT(status == SPARSE_STATUS_SUCCESS)
            << "MKLILUSolverBase failed its forward sweep at iteration "
            << i_sweep << " when calling mkl_sparse_d_trsv returned "
            << TranslateMKLSparseReturnCode(status);
        KRATOS_CATCH("")

        // Backward substitution.
        KRATOS_TRY
        lu_properties.mode = SPARSE_FILL_MODE_UPPER;
        lu_properties.diag = SPARSE_DIAG_NON_UNIT;
        const auto status = mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE,
                                              mpImpl->mRelaxation,
                                              *mpImpl->mpLUAdaptor,
                                              lu_properties,
                                              mpImpl->mBuffer.data(),
                                              SolutionView.it_begin);
        KRATOS_ERROR_IF_NOT(status == SPARSE_STATUS_SUCCESS)
            << "MKLILUSolverBase failed its backward sweep at iteration "
            << i_sweep << " when calling mkl_sparse_d_trsv returned "
            << TranslateMKLSparseReturnCode(status);
        KRATOS_CATCH("")
    }

    return true;
}


template <class TSparse, class TDense>
int MKLILUSolverBase<TSparse,TDense>::GetIterations() const noexcept
{
    return mpImpl->mIterations;
}


template <class TSparse, class TDense>
MKLILU0Solver<TSparse,TDense>::MKLILU0Solver()
    : MKLILU0Solver(Parameters())
{
}


template <class TSparse, class TDense>
MKLILU0Solver<TSparse,TDense>::MKLILU0Solver(Parameters Settings)
    : Base(Settings)
{
}


template <class TSparse, class TDense>
Parameters MKLILU0Solver<TSparse,TDense>::GetDefaultParameters()
{
    Parameters defaults(R"({
        "solver_type" : "mkl_ilu0"
    })");
    defaults.AddMissingParameters(Base::GetDefaultParameters());
    return defaults;
}


template <class TSparse, class TDense>
void MKLILU0Solver<TSparse,TDense>::Factorize(std::vector<int>& rRowExtents,
                                              std::vector<int>& rColumnIndices,
                                              std::vector<typename TSparse::DataType>& rEntries,
                                              typename Base::CSRView LhsView,
                                              const std::array<int,128>& rIntegerSettings,
                                              const std::array<typename TSparse::DataType,128>& rNumericSettings)
{
    KRATOS_TRY
    const auto ilu0_entry_count = LhsView.entry_count;

    rRowExtents.resize(LhsView.row_count + 1, 0);
    rColumnIndices.resize(ilu0_entry_count, 0.0);
    rEntries.resize(ilu0_entry_count, 0.0);

    MKL_INT return_code;
    dcsrilu0(&LhsView.row_count,
             LhsView.it_entry_begin,
             LhsView.it_row_begin,
             LhsView.it_column_begin,
             rEntries.data(),
             rIntegerSettings.data(),
             rNumericSettings.data(),
            &return_code);

    std::copy(LhsView.it_row_begin,
              LhsView.it_row_begin + LhsView.row_count + 1,
              rRowExtents.data());

    std::copy(LhsView.it_column_begin,
              LhsView.it_column_begin + LhsView.entry_count,
              rColumnIndices.data());

    KRATOS_ERROR_IF(return_code) << "dcsrilut terminated with error code " << return_code;
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
MKLILUTSolver<TSparse,TDense>::MKLILUTSolver()
    : MKLILUTSolver(Parameters())
{
}


template <class TSparse, class TDense>
MKLILUTSolver<TSparse,TDense>::MKLILUTSolver(Parameters Settings)
    : Base(Settings)
{
    KRATOS_TRY
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mFactorizationTolerance = Settings["factorization_tolerance"].Get<double>();
    mFillFactor = Settings["fill_factor"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
Parameters MKLILUTSolver<TSparse,TDense>::GetDefaultParameters()
{
    Parameters defaults(R"({
        "solver_type" : "mkl_ilut",
        "fill_factor" : 2,
        "factorization_tolerance" : 1e-2
    })");
    defaults.AddMissingParameters(Base::GetDefaultParameters());
    return defaults;
}


template <class TSparse, class TDense>
void MKLILUTSolver<TSparse,TDense>::Factorize(std::vector<int>& rRowExtents,
                                              std::vector<int>& rColumnIndices,
                                              std::vector<typename TSparse::DataType>& rEntries,
                                              typename Base::CSRView LhsView,
                                              const std::array<int,128>& rIntegerSettings,
                                              const std::array<typename TSparse::DataType,128>& rNumericSettings)
{
    KRATOS_TRY
    const auto ilut_entry_count = (2 * mFillFactor + 1) * LhsView.row_count
                                - mFillFactor * (mFillFactor + 1)
                                + 1;

    rRowExtents.resize(LhsView.row_count + 1, 0);
    rColumnIndices.resize(ilut_entry_count, 0.0);
    rEntries.resize(ilut_entry_count, 0.0);

    MKL_INT return_code;
    dcsrilut(&LhsView.row_count,
             LhsView.it_entry_begin,
             LhsView.it_row_begin,
             LhsView.it_column_begin,
             rEntries.data(),
             rRowExtents.data(),
             rColumnIndices.data(),
             &mFactorizationTolerance,
             &mFillFactor,
             rIntegerSettings.data(),
             rNumericSettings.data(),
            &return_code);

    KRATOS_ERROR_IF(return_code) << "dcsrilut terminated with error code " << return_code;
    KRATOS_CATCH("")
}


template class MKLILUSolverBase<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;


template class MKLILU0Solver<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;


template class MKLILUTSolver<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;


} // namespace Kratos
