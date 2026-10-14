//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//
// Kratos extension injected into Eigen::MatrixBase through the
// EIGEN_MATRIXBASE_PLUGIN mechanism (defined globally by the Kratos CMake
// when KRATOS_LINEAR_ALGEBRA_BACKEND is "eigen").
//
// This file is textually inserted INSIDE the Eigen::MatrixBase class body, so
// it must not have include guards, cannot include headers and must only rely
// on what Eigen itself has already made available at that point.
//
// It adds the uBLAS-style extent queries and clear() to EVERY dense Eigen
// expression (blocks, maps, products, transposes, ...), so generic Kratos code
// written against the uBLAS member surface (rA.size1(), rV.clear(), ...) also
// works on unevaluated expression arguments. The counts are returned as
// std::size_t, as uBLAS does, so the pervasive comparisons against unsigned
// counters stay warning-free.

/**
 * @brief uBLAS-style number of rows
 * @details Returns the number of rows as a std::size_t, as the uBLAS size1() accessor does.
 * @return Number of rows, as a std::size_t.
 */
inline std::size_t size1() const { return static_cast<std::size_t>(derived().rows()); }

/**
 * @brief uBLAS-style number of columns
 * @details Returns the number of columns as a std::size_t, as the uBLAS size2() accessor does.
 * @return Number of columns, as a std::size_t.
 */
inline std::size_t size2() const { return static_cast<std::size_t>(derived().cols()); }

/**
 * @brief uBLAS-style clear
 * @details Sets every entry to zero and keeps the size, as the uBLAS clear() does.
 */
inline void clear() { derived().setZero(); }

/**
 * @brief Coefficient accessors of the Eigen base
 * @details Brings Eigen's own operator()(i) and operator()(i, j) into scope next to the closure-call overloads below.
 */
using Base::operator();
/**
 * @brief uBLAS closure call on a const expression
 * @details Returns the concrete expression, as e() does on the boost::numeric::ublas::vector_expression and
 * matrix_expression bases.
 * @return Const reference to the concrete expression.
 */
inline const Derived& operator()() const { return derived(); }
/**
 * @brief uBLAS closure call
 * @details Returns the concrete expression, as e() does on the boost::numeric::ublas::vector_expression and
 * matrix_expression bases.
 * @return Reference to the concrete expression.
 */
inline Derived& operator()() { return derived(); }

/**
 * @brief Free-function form of noalias for use inside Eigen-derived classes
 * @details Called unqualified from inside a member function of a class that derives from an Eigen type (array_1d,
 * Point, Node, EigenMatrix, ...). There, class-scope lookup finds Eigen's member noalias() first and hides
 * Kratos::noalias (member lookup also suppresses ADL), so noalias(this->Coordinates()) = ... would not compile. This
 * overload accepts the argument and forwards it to Kratos::noalias. Nothing outside Eigen can be named here, so the
 * call goes through the dependent unqualified hook KratosNoAliasHook (namespace Eigen, includes/eigen_operations.h),
 * found by ADL when the call is instantiated.
 * @tparam TTarget Deduced type of the argument, a forwarding reference.
 * @param rTarget Target of the noalias call, a matrix or a matrix expression.
 * @return The noalias proxy of rTarget, as returned by Kratos::noalias.
 */
template <typename TTarget>
static inline auto noalias(TTarget&& rTarget)
{
    return KratosNoAliasHook(static_cast<TTarget&&>(rTarget));
}
