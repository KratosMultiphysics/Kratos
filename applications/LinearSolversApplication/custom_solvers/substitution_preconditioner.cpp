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

// --- External Includes ---
#include "Eigen/Sparse"

// --- Kratos LinearSolvers Includes ---
#include "custom_solvers/substitution_preconditioner.hpp"

// --- Kratos Core Includes ---
#include "spaces/ublas_space.h"


namespace Kratos {


template <class TSparse, class TDense>
struct SubstitutionPreconditioner<TSparse,TDense>::Impl {
    using EigenSparseMatrix = Eigen::SparseMatrix<
        typename TSparse::DataType,
        Eigen::RowMajor,
        std::conditional_t<
            std::is_same_v<typename TSparse::IndexType,std::size_t>,
            std::ptrdiff_t,
            void>>;

    using EigenArray = Eigen::Matrix<
        typename TSparse::DataType,
        Eigen::Dynamic,
        1>;

    using EigenSparseView = Eigen::Map<EigenSparseMatrix>;

    using EigenArrayView = Eigen::Map<EigenArray>;

    std::shared_ptr<typename TSparse::MatrixType> mpLowerTriangle, mpUpperTriangle;

    EigenSparseView mLowerTriangleView, mUpperTriangleView;

    typename TSparse::VectorType mBuffer;
}; // struct Impl


template <class TSparse, class TDense>
SubstitutionPreconditioner<TSparse,TDense>::SubstitutionPreconditioner()
    : SubstitutionPreconditioner(nullptr)
{}


template <class TSparse, class TDense>
SubstitutionPreconditioner<TSparse,TDense>::SubstitutionPreconditioner(std::shared_ptr<typename TSparse::MatrixType> pTriangle)
    : SubstitutionPreconditioner(pTriangle, nullptr)
{}


template <class TSparse, class TDense>
SubstitutionPreconditioner<TSparse,TDense>::SubstitutionPreconditioner(
    std::shared_ptr<typename TSparse::MatrixType> pLowerTriangle,
    std::shared_ptr<typename TSparse::MatrixType> pUpperTriangle)
        : mpImpl(new Impl {
            .mpLowerTriangle = pLowerTriangle,
            .mpUpperTriangle = pUpperTriangle,
            .mLowerTriangleView = typename Impl::EigenSparseView(0, 0, 0, nullptr, nullptr, nullptr),
            .mUpperTriangleView = typename Impl::EigenSparseView(0, 0, 0, nullptr, nullptr, nullptr)
        }) {
    KRATOS_TRY
        // Construct a view over the lower triangle.
        if (mpImpl->mpLowerTriangle) {
            KRATOS_ERROR_IF_NOT(mpImpl->mpLowerTriangle->index1_data().size());
            KRATOS_ERROR_IF_NOT(mpImpl->mpLowerTriangle->index2_data().size());
            KRATOS_ERROR_IF_NOT(mpImpl->mpLowerTriangle->value_data().size());
            mpImpl->mLowerTriangleView = typename Impl::EigenSparseView(
                TSparse::Size1(*mpImpl->mpLowerTriangle),
                TSparse::Size2(*mpImpl->mpLowerTriangle),
                /*only ublas space*/ mpImpl->mpLowerTriangle->value_data().size(),
                /*only ublas space*/ reinterpret_cast<std::ptrdiff_t*>(&*mpImpl->mpLowerTriangle->index1_data().begin()),
                /*only ublas space*/ reinterpret_cast<std::ptrdiff_t*>(&*mpImpl->mpLowerTriangle->index2_data().begin()),
                /*only ublas space*/ &*mpImpl->mpLowerTriangle->value_data().begin());

            mpImpl->mBuffer.resize(TSparse::Size2(*mpImpl->mpLowerTriangle));
        }

        // Construct a view over the upper triangle.
        if (mpImpl->mpUpperTriangle) {
            KRATOS_ERROR_IF_NOT(mpImpl->mpUpperTriangle->index1_data().size());
            KRATOS_ERROR_IF_NOT(mpImpl->mpUpperTriangle->index2_data().size());
            KRATOS_ERROR_IF_NOT(mpImpl->mpUpperTriangle->value_data().size());
            mpImpl->mUpperTriangleView = typename Impl::EigenSparseView(
                TSparse::Size1(*mpImpl->mpUpperTriangle),
                TSparse::Size2(*mpImpl->mpUpperTriangle),
                /*only ublas space*/ mpImpl->mpUpperTriangle->value_data().size(),
                /*only ublas space*/ reinterpret_cast<std::ptrdiff_t*>(&*mpImpl->mpUpperTriangle->index1_data().begin()),
                /*only ublas space*/ reinterpret_cast<std::ptrdiff_t*>(&*mpImpl->mpUpperTriangle->index2_data().begin()),
                /*only ublas space*/ &*mpImpl->mpUpperTriangle->value_data().begin());
        }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
SubstitutionPreconditioner<TSparse,TDense>::SubstitutionPreconditioner(SubstitutionPreconditioner&& rRhs) noexcept = default;


template <class TSparse, class TDense>
SubstitutionPreconditioner<TSparse,TDense>::~SubstitutionPreconditioner() = default;


template <class TSparse, class TDense>
SubstitutionPreconditioner<TSparse,TDense>&
SubstitutionPreconditioner<TSparse,TDense>::operator=(SubstitutionPreconditioner&& rRhs) noexcept = default;


template <class TSparse, class TDense>
void SubstitutionPreconditioner<TSparse,TDense>::Mult(
    typename TSparse::MatrixType& rLhs,
    typename TSparse::VectorType& rSolution,
    typename TSparse::VectorType& rRhs) {
        KRATOS_TRY
        TSparse::Mult(rLhs, rSolution, rRhs);
        this->ApplyLeft(rRhs);
        KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void SubstitutionPreconditioner<TSparse,TDense>::TransposeMult(
    typename TSparse::MatrixType& rLhs,
    typename TSparse::VectorType& rSolution,
    typename TSparse::VectorType& rRhs) {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(TSparse::Size(mpImpl->mBuffer) == TSparse::Size(rSolution));
        TSparse::Copy(rSolution, mpImpl->mBuffer);
        this->ApplyTransposeLeft(mpImpl->mBuffer);
        TSparse::TransposeMult(rLhs, mpImpl->mBuffer, rRhs);
        KRATOS_CATCH("")
}


template <class TSparse, class TDense>
typename TSparse::VectorType&
SubstitutionPreconditioner<TSparse,TDense>::ApplyLeft(typename TSparse::VectorType& rSolution) {
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(TSparse::Size(rSolution));
    typename Impl::EigenArrayView solution_view(
        &rSolution[0],
        TSparse::Size(rSolution));

    // Do forward substitution.
    KRATOS_ERROR_IF_NOT(mpImpl->mpLowerTriangle);
    mpImpl->mLowerTriangleView.template triangularView<Eigen::Lower>().solveInPlace(solution_view);

    // Do backward substitution.
    if (mpImpl->mpUpperTriangle) {
        // Do a backward substitution on the explicit upper triangle.
        mpImpl->mUpperTriangleView.template triangularView<Eigen::Upper>().solveInPlace(solution_view);
    } else {
        // Do a backward substitution on the lower triangle's transpose.
        mpImpl->mLowerTriangleView.template triangularView<Eigen::Upper>().solveInPlace(solution_view);
    }

    return rSolution;
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
typename TSparse::VectorType&
SubstitutionPreconditioner<TSparse,TDense>::ApplyTransposeLeft(typename TSparse::VectorType& rSolution) {
    KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented";
}


template class KRATOS_API(LINEARSOLVERS_APPLICATION) SubstitutionPreconditioner<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>>;


} // namespace Kratos
