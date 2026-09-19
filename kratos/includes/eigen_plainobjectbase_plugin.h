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
// Kratos extension injected into Eigen::PlainObjectBase through the
// EIGEN_PLAINOBJECTBASE_PLUGIN mechanism (defined globally by the Kratos
// CMake when KRATOS_LINEAR_ALGEBRA_BACKEND is "eigen").
//
// This file is textually inserted INSIDE the Eigen::PlainObjectBase class
// body, so it must not have include guards, cannot include headers and must
// only rely on what Eigen itself has already made available at that point.
//
// It adds the uBLAS-style resize overloads with an explicit "preserve" flag
// to every plain Eigen object (Eigen::Matrix, Eigen::Array, ...):
//     A.resize(rows, cols, preserve)   and   v.resize(size, preserve)
// The flag is constrained to a genuine bool so that the plain Eigen
// resize(rows, cols) never resolves to resize(size, preserve) with an
// integer-to-bool conversion.

/// uBLAS-style resize of a matrix: preserve == true keeps the existing
/// entries (Eigen's conservativeResize), preserve == false discards them.
template <typename B, typename = typename std::enable_if<std::is_same<typename std::decay<B>::type, bool>::value>::type>
inline void resize(Index rows, Index cols, B preserve)
{
    if (preserve) {
        this->conservativeResize(rows, cols);
    } else {
        this->resize(rows, cols);
    }
}

/// uBLAS-style resize of a vector: preserve == true keeps the existing
/// entries (Eigen's conservativeResize), preserve == false discards them.
template <typename B, typename = typename std::enable_if<std::is_same<typename std::decay<B>::type, bool>::value>::type>
inline void resize(Index size, B preserve)
{
    if (preserve) {
        this->conservativeResize(size);
    } else {
        this->resize(size);
    }
}
