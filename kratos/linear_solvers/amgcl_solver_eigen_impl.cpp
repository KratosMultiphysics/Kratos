//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Core includes
#include "linear_solvers/amgcl_solver.h" // AMGCLSolver
#include "spaces/eigen_space.h" // TEigenSparseSpace
#include "spaces/ublas_space.h" // TUblasDenseSpace
#include "spaces/default_spaces.h" // TDefaultDenseSpace

// See the note in amgcl_solver_impl.cpp: the implementation of AMGCLSolver is
// split between an implementation header and per-representation source files.
// This source file provides the adaptor and the instantiations for the Eigen
// sparse space.
#include "linear_solvers/amgcl_solver_impl.hpp" // AMGCLAdaptor

// STD includes
#include <optional>


namespace Kratos {


template <class TValue>
struct AMGCLAdaptor<TEigenSparseSpace<TValue>>
{
    template <int BlockSize>
    auto MakeMatrixAdaptor(const typename TEigenSparseSpace<TValue>::MatrixType& rMatrix)
    {
        if constexpr (BlockSize == 1) {
            KRATOS_TRY
            // Zero-copy view over the compressed CSR arrays of the Eigen matrix
            return amgcl::adapter::zero_copy(rMatrix.size1(),
                                             rMatrix.outerIndexPtr(),
                                             rMatrix.innerIndexPtr(),
                                             rMatrix.valuePtr());
            KRATOS_CATCH("")
        } else {
            using BlockType = amgcl::static_matrix<
                TValue,
                BlockSize,
                BlockSize
            >;

            KRATOS_TRY
            mIntermediateAdaptor.emplace(
                rMatrix.size1(),
                rMatrix.index1_data(),
                rMatrix.index2_data(),
                rMatrix.value_data());

            return amgcl::adapter::block_matrix<BlockType>(mIntermediateAdaptor.value());
            KRATOS_CATCH("")
        }
    }

    template <class TStaticMatrix>
    std::size_t BlockSystemSize(const typename TEigenSparseSpace<TValue>::MatrixType& rMatrix) const noexcept
    {
        return TEigenSparseSpace<TValue>::Size1(rMatrix) / AMGCLStaticVectorTraits<TStaticMatrix>::value;
    }

    auto MakeVectorIterator(const typename TEigenSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.size() == 0);
        return rVector.data();
    }

    auto MakeVectorIterator(typename TEigenSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.size() == 0);
        return rVector.data();
    }

private:
    using TMatrix = typename TEigenSparseSpace<TValue>::MatrixType;
    using ConstIndexProxyType = decltype(std::declval<const TMatrix&>().index1_data());
    using ConstValueProxyType = decltype(std::declval<const TMatrix&>().value_data());

    // As in the uBLAS adaptor, amgcl::adapter::block_matrix stores a reference
    // to the "matrix" (the tuple below), which must be kept alive until the
    // hierarchy construction finishes. The storage proxies are cheap views, so
    // they are stored by value.
    std::optional<std::tuple<
        std::size_t,
        ConstIndexProxyType,
        ConstIndexProxyType,
        ConstValueProxyType
    >> mIntermediateAdaptor;
}; // struct AMGCLAdaptor


// The Eigen-real AMGCLSolver instantiations exist only under the eigen
// backend: the two linear-algebra backends are mutually exclusive.
#ifdef KRATOS_USE_EIGEN_BACKEND

template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TEigenSparseSpace<double>,
    TDefaultDenseSpace<double>
>;

template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TEigenSparseSpace<float>,
    TDefaultDenseSpace<double>
>;

#endif // KRATOS_USE_EIGEN_BACKEND


} // namespace Kratos
