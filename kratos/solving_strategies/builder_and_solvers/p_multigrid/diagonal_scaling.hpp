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

#pragma once

// Project includes
#include "includes/kratos_parameters.h" // Parameters
#include "utilities/reduction_utilities.h" // CombinedReduction, MaxReduction, SumReduction

// STL includes
#include <memory> // std::unique_ptr
#include <optional> // std::optional


namespace Kratos {


/** @brief Class handling expressions involving the diagonal of a sparse matrix.
 *  @details
 *  @subsection Configuration
 *              This class can either be configured with a numeric literal that does not
 *              depend on the matrix's diagonal, or a string that is a function of:
 *              - @p "max" The infinitiy norm (abs max) of the matrix' diagonal.
 *              - @p "norm" 2-norm of the matrix' diagonal.
 *
 *  @subsection Usage
 *              Invoke @ref Scaling::Cache "Cache" with a matrix to evaluate the stored
 *              expression for later use. The computed value can be retrived from
 *              @ref Scaling::Evaluate "Evaluate". You should call
 *              @ref Scaling::ClearCache "ClearCache" whenever the input matrix changes.
 *
 *  @note The input expression is parsed and evaluated with @p tinyexpr.
 */
class KRATOS_API(KRATOS_CORE) Scaling final
{
public:
    /// @brief Construct using the default settings.
    Scaling();

    /// @brief Construct from @ref Parameters representing a numeric literal or a string.
    /// @details
    Scaling(Parameters Settings);

    /// @internal
    ~Scaling();

    /// @brief Compute and store the provided expression using the diagonal of an input matrix.
    /// @note Invoking this function is mandatory before calling @ref Scaling::Evaluate "Evaluate".
    template <class TSparse>
    void Cache(const typename TSparse::MatrixType& rMatrix);

    /// @brief Reset the expression's cached value.
    /// @note Subsequent calls to @ref Scaling::Evaluate will fail unless @ref Scaling::Cache is called beforehand.
    void ClearCache();

    /// @brief Get the cached value of the stored expression.
    /// @warning A call to @ref Scaling::Cache must precede any call to @p Evaluate.
    double Evaluate() const;

private:
    void Cache(double DiagonalMax, double DiagonalNorm);

    Scaling(Scaling&&) = delete;

    Scaling(const Scaling&) = delete;

    Scaling& operator=(Scaling&&) = delete;

    Scaling& operator=(const Scaling&) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class Scaling


/// @internal
template <class TSparse>
void NormalizeRows(typename TSparse::MatrixType& rLhs,
                   typename TSparse::VectorType& rRhs);


/// @internal
template <class TSparse>
void NormalizeSystem(typename TSparse::MatrixType& rLhs,
                     typename TSparse::VectorType& rRhs,
                     typename TSparse::DataType Coefficient);


namespace detail {


/// @internal
template <class TSparse>
std::optional<typename TSparse::DataType> FindDiagonal(const typename TSparse::MatrixType& rMatrix,
                                                       const std::size_t iRow) noexcept
{
    const auto i_entry_begin = rMatrix.index1_data()[iRow];
    const auto i_entry_end   = rMatrix.index1_data()[iRow + 1];

    const auto it_column_begin = rMatrix.index2_data().begin() + i_entry_begin;
    const auto it_column_end = rMatrix.index2_data().begin() + i_entry_end;

    // Look for the diagonal entry in the current row.
    const auto it_column = std::lower_bound(it_column_begin, it_column_end, iRow);

    return (it_column == it_column_end || *it_column != iRow)
         ? std::optional<typename TSparse::DataType>()
         : rMatrix.value_data()[rMatrix.index1_data()[iRow] + std::distance(it_column_begin, it_column)];
}


} // namespace detail


/// @internal
template <class TSparse>
void Scaling::Cache(const typename TSparse::MatrixType& rMatrix)
{
    KRATOS_TRY

    using Value = typename TSparse::DataType;
    using Reduction = CombinedReduction<MaxReduction<Value>,SumReduction<Value>>;

    Value norm_coefficient = rMatrix.size1()
        ? std::pow(static_cast<Value>(1) / static_cast<Value>(rMatrix.size1()), static_cast<Value>(2))
        : static_cast<Value>(1);

    const auto [abs_max, square_norm] = IndexPartition<typename TSparse::IndexType>(rMatrix.size1()).template for_each<Reduction>([&rMatrix, norm_coefficient](auto i_row){
        const auto maybe_diagonal_entry = detail::FindDiagonal<TSparse>(rMatrix, i_row);
        const Value diagonal_entry = maybe_diagonal_entry.has_value() ? *maybe_diagonal_entry : static_cast<Value>(0);
        return std::make_tuple(std::abs(diagonal_entry), norm_coefficient * diagonal_entry * diagonal_entry);
    });

    this->Cache(abs_max, std::sqrt(square_norm));

    KRATOS_CATCH("")
}


} // namespace Kratos
