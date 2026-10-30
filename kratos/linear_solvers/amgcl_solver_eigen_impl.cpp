//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// The Eigen adaptor and the Eigen-space AMGCLSolver instantiations exist only
// under the eigen backend (amgcl_solver_impl.cpp provides the uBLAS ones).
#ifdef KRATOS_USE_EIGEN_BACKEND

// Core includes
#include "linear_solvers/amgcl_solver.h" // AMGCLSolver
#include "spaces/default_spaces.h" // TEigenSparseSpace, TDefaultDenseSpace

// See the note in amgcl_solver_impl.cpp: the implementation of AMGCLSolver is
// split between an implementation header and per-representation source files.
// This source file provides the adaptor and the instantiations for the Eigen
// sparse space.
#include "linear_solvers/amgcl_solver_impl.hpp" // AMGCLAdaptor

// STD includes
#include <optional>


namespace Kratos {


/**
 * @brief Adaptor of the Eigen sparse space for AMGCLSolver
 * @details Specialization of AMGCLAdaptor for TEigenSparseSpace<TValue>. It converts the Eigen sparse matrices and
 * vectors into the amgcl types that AMGCLSolver works with.
 * @tparam TValue Scalar type of the Eigen sparse space.
 */
template <class TValue>
struct AMGCLAdaptor<TEigenSparseSpace<TValue>>
{
    /**
     * @brief Builds the amgcl matrix adaptor of an Eigen sparse matrix
     * @details For BlockSize == 1 it returns a zero-copy amgcl view over the compressed CSR arrays when the index type
     * has the width of std::ptrdiff_t. Otherwise it returns a tuple of the storage proxies, which amgcl widens while
     * copying. For BlockSize > 1 the size and the storage proxies are kept in mIntermediateAdaptor, and an amgcl block
     * matrix view over them is returned.
     * @tparam BlockSize Size of the blocks of the amgcl matrix, 1 for the scalar CSR view.
     * @param rMatrix Eigen sparse matrix to adapt.
     * @return amgcl adapter of the matrix: a zero-copy CSR view, a tuple of storage proxies, or a block matrix view.
     */
    template <int BlockSize>
    auto MakeMatrixAdaptor(const typename TEigenSparseSpace<TValue>::MatrixType& rMatrix)
    {
        if constexpr (BlockSize == 1) {
            KRATOS_TRY
            using StorageIndexType = std::remove_pointer_t<decltype(rMatrix.outerIndexPtr())>;
            if constexpr (sizeof(StorageIndexType) == sizeof(std::ptrdiff_t)) {
                // Zero-copy view over the compressed CSR arrays of the Eigen matrix
                return amgcl::adapter::zero_copy(rMatrix.size1(),
                                                 rMatrix.outerIndexPtr(),
                                                 rMatrix.innerIndexPtr(),
                                                 rMatrix.valuePtr());
            } else {
                // amgcl's zero_copy reinterprets the index arrays as ptrdiff_t
                // and cannot accept the (default) 32-bit indices; hand amgcl
                // the storage proxies through its generic range-tuple
                // interface instead, which widens the indices while copying
                // into the hierarchy it builds anyway.
                return std::make_tuple(rMatrix.size1(),
                                       rMatrix.index1_data(),
                                       rMatrix.index2_data(),
                                       rMatrix.value_data());
            }
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

    /**
     * @brief Number of block rows of the system
     * @details Returns the number of rows of rMatrix divided by the block size of TStaticMatrix, given by
     * AMGCLStaticVectorTraits<TStaticMatrix>::value.
     * @tparam TStaticMatrix Static matrix type whose block size is used.
     * @param rMatrix Eigen sparse matrix of the system.
     * @return Number of block rows of the system.
     */
    template <class TStaticMatrix>
    std::size_t BlockSystemSize(const typename TEigenSparseSpace<TValue>::MatrixType& rMatrix) const noexcept
    {
        return TEigenSparseSpace<TValue>::Size1(rMatrix) / AMGCLStaticVectorTraits<TStaticMatrix>::value;
    }

    /**
     * @brief Read-only iterator over the entries of a vector
     * @details Raises an error for an empty vector. Returns a const pointer to the first entry of the vector storage.
     * @param rVector Vector whose entries are visited.
     * @return Read-only pointer to the first entry of the vector storage.
     */
    auto MakeVectorIterator(const typename TEigenSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.size() == 0);
        return rVector.data().begin();
    }

    /**
     * @brief Iterator over the entries of a vector
     * @details Raises an error for an empty vector. Returns a pointer to the first entry of the vector storage, so the
     * entries can be written through it.
     * @param rVector Vector whose entries are visited.
     * @return Pointer to the first entry of the vector storage.
     */
    auto MakeVectorIterator(typename TEigenSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.size() == 0);
        return rVector.data().begin();
    }

private:
    /// Eigen sparse matrix type of the space. Alias of TEigenSparseSpace<TValue>::MatrixType.
    using TMatrix = typename TEigenSparseSpace<TValue>::MatrixType;
    /// Index array view of a const matrix. Type returned by index1_data() on a const matrix. The same type is used for
    /// both index arrays.
    using ConstIndexProxyType = decltype(std::declval<const TMatrix&>().index1_data());
    /// Value array view of a const matrix. Type returned by value_data() on a const matrix.
    using ConstValueProxyType = decltype(std::declval<const TMatrix&>().value_data());

    // As in the uBLAS adaptor, amgcl::adapter::block_matrix stores a reference
    // to the "matrix" (the tuple below), which must be kept alive until the
    // hierarchy construction finishes. The storage proxies are cheap views, so
    // they are stored by value.
    /**
     * @brief Intermediate tuple of the block matrix adaptor
     * @details Holds the matrix size and the storage proxies that amgcl::adapter::block_matrix refers to, so that they
     * stay alive until the hierarchy is built. It is only filled in the block-size branch of MakeMatrixAdaptor.
     */
    std::optional<std::tuple<
        std::size_t,
        ConstIndexProxyType,
        ConstIndexProxyType,
        ConstValueProxyType
    >> mIntermediateAdaptor;
}; // struct AMGCLAdaptor


template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TEigenSparseSpace<double>,
    TDefaultDenseSpace<double>
>;

template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TEigenSparseSpace<float>,
    TDefaultDenseSpace<double>
>;


} // namespace Kratos

#endif // KRATOS_USE_EIGEN_BACKEND
