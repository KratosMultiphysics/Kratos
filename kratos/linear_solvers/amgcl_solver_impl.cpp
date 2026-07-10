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
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace

// The implementation of AMGCLSolver is split between
// - an implementation header ("linear_solvers/amgcl_solver_impl.hpp")
// - implementation sources (e.g.: this source file)
//
// The reason is twofold:
// - includes from the AMGCL library are extremely heavy, so they are
//   avoided in the class declaration ("linear_solvers/amgcl_solver.h").
//   Instead, the implementation header includes them and defines logic
//   common to any matrix/vector representations. Each source file that
//   defines an instantiation of AMGCLSolver includes the implementation
//   header.
// - Shared memory and distributed memory matrix/vector representations
//   are handled in separate source files to avoid adding a Trilinos
//   dependency to core.
#include "linear_solvers/amgcl_solver_impl.hpp" // AMGCLAdaptor

// STD includes
#include <optional>


namespace Kratos {


template <class TValue>
struct AMGCLAdaptor<TUblasSparseSpace<TValue>>
{
    template <int BlockSize>
    auto MakeMatrixAdaptor(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix)
    {
        if constexpr (BlockSize == 1) {
            KRATOS_TRY
            return amgcl::adapter::zero_copy(rMatrix.size1(),
                                             rMatrix.index1_data().begin(),
                                             rMatrix.index2_data().begin(),
                                             rMatrix.value_data().begin());
            KRATOS_CATCH("")
        } else {
            using BlockType = amgcl::static_matrix<
                TValue,
                BlockSize,
                BlockSize
            >;

            KRATOS_TRY
            mIntermediateAdaptor.emplace(std::tuple_cat(
                std::tuple<std::size_t>(rMatrix.size1()),
                std::tie(rMatrix.index1_data(),
                         rMatrix.index2_data(),
                         rMatrix.value_data())
            ));

            return amgcl::adapter::block_matrix<BlockType>(mIntermediateAdaptor.value());
            KRATOS_CATCH("")
        }
    }

    template <class TStaticMatrix>
    std::size_t BlockSystemSize(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix) const noexcept
    {
        return TUblasSparseSpace<TValue>::Size1(rMatrix) / AMGCLStaticVectorTraits<TStaticMatrix>::value;
    }

    auto MakeVectorIterator(const typename TUblasSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.empty());
        return &*rVector.begin();
    }

    auto MakeVectorIterator(typename TUblasSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.empty());
        return &*rVector.begin();
    }

private:
    using TMatrix = typename TUblasSparseSpace<TValue>::MatrixType;

    // amgcl::adapter::block_matrix constructs a class
    // that stores a reference to the "matrix" passed
    // into it, which in this case means the tuple
    // defined below. We need to keep it alive until
    // the hierarchy construction finishes, hence the
    // convoluted member variable.
    // Optional is used here to represent an invalid state
    // of the matrix view, before InitializeSolutionStep is
    // called.
    std::optional<std::tuple<
        std::size_t,
        const typename TMatrix::index_array_type&,
        const typename TMatrix::index_array_type&,
        const typename TMatrix::value_array_type&
    >> mIntermediateAdaptor;
}; // struct AMGCLAdaptor


template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>
>;

template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos
