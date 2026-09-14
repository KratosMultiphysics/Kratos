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

/// uBLAS-style number of rows.
inline std::size_t size1() const { return static_cast<std::size_t>(derived().rows()); }

/// uBLAS-style number of columns.
inline std::size_t size2() const { return static_cast<std::size_t>(derived().cols()); }

/// uBLAS-style clear(): zero every entry, keeping the size.
inline void clear() { derived().setZero(); }

/// uBLAS closure call: e() yields the concrete expression, as on the
/// boost::numeric::ublas::vector_expression / matrix_expression bases. The
/// using-declaration keeps Eigen's own coefficient accessors operator()(i)
/// and operator()(i, j) visible next to the added overload.
using Base::operator();
inline const Derived& operator()() const { return derived(); }
inline Derived& operator()() { return derived(); }
