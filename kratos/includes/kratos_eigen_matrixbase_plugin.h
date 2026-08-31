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
//
// Kratos extension injected into Eigen::MatrixBase through the
// EIGEN_MATRIXBASE_PLUGIN mechanism (defined globally by the Kratos CMake
// when KRATOS_LINEAR_ALGEBRA_BACKEND is "eigen").
//
// It adds implicit conversion of ANY dense Eigen expression (lazy products,
// sums, transposes, blocks, ...) to the dynamically-sized
// boost::numeric::ublas dense containers, restoring the uBLAS-backend
// semantics of idioms such as
//     Vector v = prod(rRotation, rLocalVector);
//     Matrix m = trans(rRotation);
//     rUblasVector = rArray + rOtherArray;
// where the target remains a uBLAS type under both backends.
//
// This file is textually inserted INSIDE the Eigen::MatrixBase class body, so
// it cannot include headers and must not name the boost namespace (Eigen may
// be included before boost in a translation unit). The uBLAS targets are
// therefore identified STRUCTURALLY: only the boost dense containers have the
// const_closure_type / array_type-with-allocator nested types required below,
// so no other conversion is affected.

/// uBLAS-style extent queries on EVERY dense Eigen expression (blocks, maps,
/// products, transposes, ...), not only on the concrete Kratos wrapper types
/// (whose identical std::size_t-returning members shadow these). Adopted from
/// external review feedback: it removes a whole class of per-site shims where
/// generic code calls size1()/size2() on an unevaluated expression argument.
inline std::size_t size1() const { return static_cast<std::size_t>(derived().rows()); }
inline std::size_t size2() const { return static_cast<std::size_t>(derived().cols()); }

/// uBLAS-style clear() (value zeroing) for writable dense expressions.
inline void clear() { derived().setZero(); }

/// Implicit conversion of a (compile-time) vector-shaped expression to a
/// dynamically-allocated uBLAS dense vector (Kratos::Vector and friends).
template<class TUblasDenseVector>
operator TUblasDenseVector() const
    requires ((Eigen::internal::traits<Derived>::RowsAtCompileTime == 1 ||
               Eigen::internal::traits<Derived>::ColsAtCompileTime == 1) &&
              requires(TUblasDenseVector target, std::size_t index) {
                  typename TUblasDenseVector::const_closure_type;             // a boost::numeric::ublas container...
                  typename TUblasDenseVector::array_type::allocator_type;     // ...with dynamically allocated storage
                  target.resize(index, false);
                  { target[index] };
              } &&
              !requires(TUblasDenseVector target) { target.size1(); } &&      // not a matrix
              std::is_convertible_v<Scalar, typename TUblasDenseVector::value_type>)
{
    const auto evaluated = derived().eval();
    using target_size_type = typename TUblasDenseVector::size_type;
    TUblasDenseVector result(static_cast<target_size_type>(evaluated.size()));
    for (Index i = 0; i < evaluated.size(); ++i) {
        result[static_cast<target_size_type>(i)] = evaluated(i);
    }
    return result;
}

/// Implicit conversion of a matrix-shaped expression to a dynamically
/// allocated uBLAS dense matrix (Kratos::Matrix and friends).
template<class TUblasDenseMatrix>
operator TUblasDenseMatrix() const
    requires (Eigen::internal::traits<Derived>::RowsAtCompileTime != 1 &&
              Eigen::internal::traits<Derived>::ColsAtCompileTime != 1 &&
              requires(TUblasDenseMatrix target, std::size_t index) {
                  typename TUblasDenseMatrix::const_closure_type;             // a boost::numeric::ublas container...
                  typename TUblasDenseMatrix::array_type::allocator_type;     // ...with dynamically allocated storage
                  typename TUblasDenseMatrix::orientation_category;           // ...that is a matrix
                  { target(index, index) };
              } &&
              std::is_convertible_v<Scalar, typename TUblasDenseMatrix::value_type>)
{
    const auto evaluated = derived().eval();
    using target_size_type = typename TUblasDenseMatrix::size_type;
    TUblasDenseMatrix result(static_cast<target_size_type>(evaluated.rows()),
                             static_cast<target_size_type>(evaluated.cols()));
    for (Index i = 0; i < evaluated.rows(); ++i) {
        for (Index j = 0; j < evaluated.cols(); ++j) {
            result(static_cast<target_size_type>(i), static_cast<target_size_type>(j)) = evaluated(i, j);
        }
    }
    return result;
}
