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
#include "tinyexpr/tinyexpr.h" // te_expr, te_compile, te_eval, te_free

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // DiagonalScaling, ParseDiagonalScaling, GetDiagonalScaleFactor
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "utilities/reduction_utilities.h" // MaxReduction, AbsMaxReduction
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE

// STL includes
#include <optional> // std::optional


namespace Kratos {


struct Scaling::Impl {
    using ExpressionPtr = std::unique_ptr<te_expr,std::function<void(te_expr*)>>;

    std::optional<double> mMaybeAbsMax;

    std::optional<double> mMaybeNorm;

    mutable ExpressionPtr mpExpression;
}; // struct Scaling::Impl


Scaling::Scaling()
    : Scaling(Parameters())
{
}


Scaling::Scaling(Parameters Settings)
    : mpImpl(new Impl {/*mMaybeAbsMax=*/0.0,
                       /*mMaybeNorm=*/  0.0,
                       /*mpExpression*/ Impl::ExpressionPtr(nullptr, [](te_expr*){})})
{
    // Sanity checks.
    KRATOS_ERROR_IF_NOT(Settings.IsNumber() || Settings.Is<std::string>())
        << "Expecting a numeric literal or a string, but got " << Settings;

    KRATOS_TRY

    std::string expression;
    if (Settings.IsNumber()){
        expression = std::to_string(Settings.Get<double>());
    } /*if Settings.IsNumber*/ else {
        expression = Settings.Get<std::string>();
    } /*if Settings.Is<std::string>*/

    te_variable variables[] = {
        {"max",  &mpImpl->mMaybeAbsMax.value()},
        {"norm", &mpImpl->mMaybeNorm.value()}
    };

    int error = 0;
    mpImpl->mpExpression = Impl::ExpressionPtr(/*function ptr: */te_compile(expression.c_str(), variables, 2, &error),
                                               /*deleter: */     [](te_expr* p){if (p) te_free(p);});
    KRATOS_ERROR_IF_NOT(mpImpl->mpExpression.get())
        << "Failed to compile expression (error at position " << error - 1 << " in \"" << expression << "\"). "
        << "Make sure the expression is compatible with TinyExpr and is a function of "
        << "only \"max\" and \"norm\".";

    this->ClearCache();
    KRATOS_CATCH("")
}


Scaling::~Scaling() = default;


void Scaling::ClearCache()
{
    mpImpl->mMaybeAbsMax.reset();
    mpImpl->mMaybeNorm.reset();
}


double Scaling::Evaluate() const
{
    KRATOS_ERROR_IF_NOT(mpImpl->mMaybeAbsMax.has_value() && mpImpl->mMaybeNorm.has_value())
        << "Matrix diagonal norms are not cached. Make sure to invoke Scaling::Cache before Scaling::Evaluate.";

    return te_eval(mpImpl->mpExpression.get());
}


void Scaling::Cache(double DiagonalMax, double DiagonalNorm)
{
    mpImpl->mMaybeAbsMax = DiagonalMax;
    mpImpl->mMaybeNorm = DiagonalNorm;
}


template <class TSparse>
void NormalizeRows(typename TSparse::MatrixType& rLhs,
                   typename TSparse::VectorType& rRhs)
{
    KRATOS_TRY
    IndexPartition<std::size_t>(rLhs.size1()).for_each([&rLhs, &rRhs](const std::size_t i_row){
        const auto maybe_diagonal_entry = detail::FindDiagonal<TSparse>(rLhs, i_row);
        KRATOS_ERROR_IF_NOT(maybe_diagonal_entry.has_value()) << "row " << i_row << " has no diagonal entry";
        const auto diagonal_entry = maybe_diagonal_entry.value();

        if (diagonal_entry) {
            const auto scale = static_cast<typename TSparse::DataType>(1) / diagonal_entry;

            const auto i_entry_begin = rLhs.index1_data()[i_row];
            const auto i_entry_end = rLhs.index1_data()[i_row + 1];
            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                rLhs.value_data()[i_entry] *= scale;
            } // for i_entry in range(i_entry_begin, i_entry_end)

            rRhs[i_row] *= scale;
        } /*if diagonal_entry*/ else {
            const auto i_entry_begin = rLhs.index1_data()[i_row];
            const auto i_entry_end = rLhs.index1_data()[i_row + 1];
            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = rLhs.index2_data()[i_entry];
                const auto value = rLhs.value_data()[i_entry];
                KRATOS_ERROR_IF(i_column != i_row && value) << "the diagonal of row " << i_row << " vanishes, but has off-diagonal components";
                //if (i_column == i_row) rLhs.value_data()[i_entry] = static_cast<typename TSparse::DataType>(1);
            } // for i_entry in range(i_entry_begin, i_entry_end)
            //KRATOS_ERROR_IF(rRhs[i_row]) << "row " << i_row << " is empty, but the corresponding RHS component is " << rRhs[i_row];
        } // if not diagonal_entry
    });
    KRATOS_CATCH("")
}


template <class TSparse>
void NormalizeSystem(typename TSparse::MatrixType& rLhs,
                     typename TSparse::VectorType& rRhs,
                     typename TSparse::DataType Coefficient)
{
    block_for_each(&*rLhs.value_data().begin(),
                   (&*rLhs.value_data().begin()) + rLhs.value_data().size(),
                  [Coefficient](typename TSparse::DataType& r_entry){r_entry *= Coefficient;});
    block_for_each(&*rRhs.begin(),
                   (&*rRhs.begin()) + rRhs.size(),
                  [Coefficient](typename TSparse::DataType& r_entry){r_entry *= Coefficient;});
}


// Explicit instantiations for UBLAS.
#define KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS(TSparseSpace)                                                        \
    template void NormalizeRows<TSparseSpace>(typename TSparseSpace::MatrixType&,                                   \
                                              typename TSparseSpace::VectorType&);                                  \
    template void NormalizeSystem<TSparseSpace>(typename TSparseSpace::MatrixType&,                                 \
                                                typename TSparseSpace::VectorType&,                                 \
                                                typename TSparseSpace::DataType)


KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS(TUblasSparseSpace<double>);

KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS(TUblasSparseSpace<float>);

#undef KRATOS_DEFINE_DIAGONAL_SCALE_FUNCTIONS


} // namespace Kratos
